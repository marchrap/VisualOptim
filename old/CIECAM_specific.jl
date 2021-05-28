# Necessary imports
using BenchmarkTools
using ForwardDiff
using StaticArrays
using LinearAlgebra
using ProgressMeter

@inline function rgb2linear(V)
    #=
        An inline function that converts the RGB values to their linear counterparts

        Params:
            V: the original RGB vector.
    =#
    v = V / 255.

    if v <= 0.04045
        v *= 25.0 / 323.0
    else
        v = ((200.0 * v + 11.0) / 211.0)^2.4
    end

    return v
end

function CAM16Model()
    #=
        Function that generates the function resposible for converting a point to
        the CIECAM16 UCS space.
    =#

    # Define model parameters
    X_w = [95.047; 100.0; 108.883]
    Y_b = 20.0
    X_wr = Y_wr = Z_wr = 100.0
    L_A = 40.0
    F = 1.0
    c = 0.69
    N_c = 1.0
    M_16 = [0.401288 0.650173 -0.051461;
            -0.250268 1.204414 0.045854;
            -0.002079 0.048952 0.953127]

    # Calculate all values independent of the input
    R_w = M_16 * X_w
    D = F * (1 - 1 / 3.6 * exp((-L_A - 42) / 92))
    if D > 1
        D = 1
    elseif D < 0
        D = 0
    end
    D_R = D * X_w[2] ./ R_w .+ 1.0 .- D
    k = 1 / (5 * L_A + 1)
    F_L = k^4 * L_A + 0.1 * (1 - k^4)^2 * (5 * L_A)^(1 / 3)
    n = Y_b / X_w[2]
    z_w = 1.48 + sqrt(n)
    N_bb = 0.725 / (n^0.2)
    N_cb = N_bb
    R_wc = D_R .* R_w
    R_aw = 400.0 .* ((F_L .* R_wc) ./ 100.0).^0.42 ./ (((F_L .* R_wc) ./ 100.0).^0.42 .+ 27.13)
    A_w = (2 * R_aw[1] + R_aw[2] + 1 / 20 * R_aw[3]) * N_bb

    let D_G = D_R[2], D_B = D_R[3], D_R = D_R[1], n = n, F_L = F_L, z_w = z_w,
          c = c, N_bb = N_bb, A_w = A_w

        return  function rgb2cam(V_r, V_g, V_b)
                    #=
                        Converts the RGB values to CIECAM16 UCS values.

                        Params:
                            V_r: The red component of the color.
                            V_g: The green component of the color.
                            V_b: The blue component of the color.
                    =#
                    # Move RGB values to linear ones
                    r = rgb2linear(V_r)
                    g = rgb2linear(V_g)
                    b = rgb2linear(V_b)

                    # Convert linear RBG values to XYZ
                    x = 0.4124564 * r + 0.3575761 * g + 0.1804375 * b
                    y = 0.2126729 * r + 0.7151522 * g + 0.0721750 * b
                    z = 0.0193339 * r + 0.1191920 * g + 0.9503041 * b

                    # TODO change namign to be correct with the paper
                    # Calculate cone responses
                    R_y = 0.401288 * x + 0.650173 * y + -0.051461 * z
                    G_y = -0.250268 * x + 1.204414 * y + 0.045854 * z
                    B_y = -0.002079 * x + 0.048952 * y + 0.953127 * z

                    # Color adaptation
                    R_c = D_R * R_y
                    G_c = D_G * G_y
                    B_c = D_B * B_y

                    # Postadaptation
                    R_a = 400.0 * sign(R_c) * (F_L * abs(R_c))^0.42 / ((F_L * abs(R_c))^0.42 + 27.13)
                    G_a = 400.0 * sign(G_c) * (F_L * abs(G_c))^0.42 / ((F_L * abs(G_c))^0.42 + 27.13)
                    B_a = 400.0 * sign(B_c) * (F_L * abs(B_c))^0.42 / ((F_L * abs(B_c))^0.42 + 27.13)

                    # Components
                    a = R_a - 12.0 / 11.0 * G_a + 1.0 / 11.0 * B_a
                    b = 1.0 / 9.0 * R_a + 1.0 / 9.0 * G_a - 2.0 / 9.0 * B_a
                    h = atan(b, a)

                    # Eccentricity
                    e_t = 1.0 / 4.0 * (cos(h + 2.0) + 3.8)

                    # Chromatic resposne
                    A = (2.0 * R_a + G_a + 1.0 / 20.0 * B_a) * N_bb

                    # Lightness
                    J = 100.0 * (A / A_w)^(c * z_w)

                    # Colorfullness
                    t = 50000.0 / 13.0 * N_c * N_bb * e_t * sqrt(a^2 + b^2) / (R_a + G_a + 21.0 / 20.0 * B_a + 0.305)
                    α = t^0.9 * (1.64 - 0.29^n)^0.73
                    C1 = α * sqrt(J / 100.0)
                    M = C1 * F_L^0.25

                    # Move to UCS
                    J_prime = 1.7 * J / (1 + 0.007 * J)
                    M_prime = log(1 + 0.0228 * M) / 0.0228
                    a_prime = M_prime * cos(h)
                    b_prime = M_prime * sin(h)

                    # Return the result
                    return J_prime, a_prime, b_prime
        end
    end
end

function differentiator(distance, config)
    #=
        The function that generates a function that evaluates the hessian, gradient
        and the function value.

        Params:
            distance: The distance metric.
            config: The ForwardDiff config of the hessian calculation
    =#
    # Make the closure more performant
    let distance = distance, config = config
        # Return the hessian of the function according to the distance of the
        # point to the given v. Here we are using DiffResult for performance.
        return (u, result) -> ForwardDiff.hessian!(result, distance, u, config)
    end
end

function distance(V_r::Float64, V_g::Float64, V_b::Float64, eps::Float64,
     rgb2cam::Function)
    #=
        Generates a distance function from the provided point that can be reused.

        Params:
            V_r: The red part of the color.
            V_g: The green part of the color.
            V_b: The blue part of the color.
            Eps: The maximum allowable deviation from the original point.
    =#
    # Obtain the point v coordinates in the space
    v_j, v_a, v_b = rgb2cam(V_r, V_g, V_b)

    # Make the closure performant
    let v_j = v_j, v_a = v_a, v_b = v_b, eps = eps
        return  function et(U)
                    #=
                        Evaluates the distance between an RGB point and the
                        original point used to generate this function.

                        Params:
                            U: The RGB triplet used to evaluate distance to the
                               point used to generate this function.
                    =#
                    # TODO add the boundary function for the positivity contraints
                    # Check if any point is out of the possible RGB range
                    # and reject the point if so
                    if U[1] < 0 || U[2] < 0 || U[3] < 0
                        return Inf
                    end

                    # Obtain the point position of the color in the CIECAM16
                    u_j, u_a, u_b = rgb2cam(U[1], U[2], U[3])

                    # Calculate ΔE
                    delta_E = ((v_j - u_j)^2 + (v_a - u_a)^2 + (v_b - u_b)^2)
                    delta = 1.41 * delta_E^(0.63 * 0.5)
                    # If ΔE is zero then the square root's reivative is undefinded
                    # And we should not use the function but return just the original
                    # function
                    #println(eps - delta)
                    # return ((U[1]/255)^2 + (U[2]/255)^2 + (U[3]/255)^2) + max(0, delta-eps+1.0)^100
                    if delta_E == 0
                        return 10^2 * (U[1]^2 + U[2]^2 + U[3]^2) / (255.0^2) - log(eps)

                    # If the constraint is broken we output infinity
                    elseif eps < delta
                        return Inf

                    # If everything is correct we output the exact function required
                    else
                        return 10^2 * (U[1]^2 + U[2]^2 + U[3]^2) / (255.0^2) - log(eps - delta) #- log(U[1]) - log(U[2]) - log(U[3])
                    end
        end
    end
end

function backtracking!(x, delta_x, f_value, grad_value, f)
    #=
        The implementation of the backtracking algorithm.

        Params:
            x: The x point.
            delta_x: The delta_x direction. Will be mutated.
            f_value: The value of the funciton at the current point.
            grad_value: The value of the gradient dotted with delta_x.
    =#
    # Define the initial parameters
    alpha = 0.2
    beta = 0.5
    t = 1.0

    # Initialize the Δx with the initial t
    delta_x .*= t

    # Define number of iterations
    iters = 0

    println("F: ", f(x + delta_x))
    println("right: ", f_value + alpha * t * grad_value)

    # Run standard backtracking
    while f(x + delta_x) > f_value + alpha * t * grad_value && iters < 1000
        # Update t
        t *= beta

        # Update Δx
        delta_x .*= beta

        # println("F: ", f(x + delta_x))
        # println("right: ", f_value + alpha * t * grad_value)
        # println("value: ", x, " delta: ", delta_x, " total: ", x + delta_x)

        # Increase number of iterations
        iters += 1
    end
end

function Newton!(x, result, differentiate!, f, tol)
    #=
        The function that runs the Newton algorithm.

        Params:
            x: The starting point. Will be mutated.
            result: The DiffResult that will be used to store the results of
                    automatic differetiation.
            f: The function that will be optimized.
            tol: The tolerance to which the algorithm should run.
    =#
    # Define the initial Δx
    delta_x = zeros(MVector{3})

    # Define number of iterations
    iters = 0

    # Run the loop
    while iters < 1000
        # Calculate function value, gradient and hessian
        result = differentiate!(x, result)

        # Obtain the stopping epsilon and the delta_x
        delta_x .= - DiffResults.hessian(result) \ DiffResults.gradient(result)
        eps = dot(DiffResults.gradient(result), delta_x)
        println(x)
        println(delta_x)
        println(DiffResults.gradient(result))
        println(DiffResults.hessian(result))
        println("diff: ", DiffResults.value(result))
        println("f direct: ", f(x))
        println("constant: ", f([0.38125, 0.38125, 0.006249999999999997]))

        # Deal with the situation when the gradient or hessian go to a NaN
        if isnan(eps)
            break
        end

        # Calculate the step and apply it to current x
        backtracking!(x, delta_x, DiffResults.value(result), eps, f)
        x .+= delta_x
        println(delta_x)
        println(eps^2/2)
        println(tol)
        println(eps^2 / 2 < tol)
        println(abs.(delta_x))
        println(any(abs.(delta_x) .<= 1e-10))
        println(x)

        # Check the stopping condition
        if eps^2 / 2 < tol || all(abs.(delta_x) .<= 1e-10)
            break
        end

        # Increase the number of iterations
        iters += 1
    end
end

function optimize!(V, rgb2cam, config, result)
    #=
        The optimization function that is run to obtain the optimal results.

        Params:
            V: The original vector of RGB that will be mutated to the optimal one.
            rgb2cam: The conversion function that will be used to move the points
                     to the RGB to CIECAM16 space.
            config: The config of the ForwardDiff framework.
            result: The result variable to store the differentiation results.
    =#
    # Get the distance function
    dist = distance(V[1], V[2], V[3], 2.0, rgb2cam)

    # Obtain the differentiation function
    differentiate! = differentiator(dist, config)

    # Run the Newton's algorithm
    Newton!(V .+ 0.1, result, differentiate!, dist, 1e-10)
end

function optimizeSpace!(space)
    #=
        Used to optimize a whole array of values.

        Params:
            space: The array of RGB triples that will be optimized.
    =#
    # Define model
    rgb2cam = CAM16Model()

    # Define an example point
    t = [150.0, 130.0, 113.0]

    # Define the distance and example configurations that can be reused
    dist = distance(t[1], t[2], t[3], 1.0, rgb2cam)
    result = DiffResults.HessianResult(t)
    config = ForwardDiff.HessianConfig(dist, result, t, ForwardDiff.Chunk{3}())

    N = length(A)
    p = Progress(N, 1, "Computing the matrix... ", 50)
    update!(p, 0)
    jj = Threads.Atomic{Int}(0)
    l = Threads.SpinLock()

    # Run over the points
    for i in space
        print(i)
        dist = distance(i[1], i[2], i[3], 2.0, rgb2cam)
        solution = Optim.minimizer(optimize(dist, i))
        i .= solution #optimize!(i, rgb2cam, config, result)
        Threads.atomic_add!(jj, 1)
        Threads.lock(l)
        update!(p, jj[])
        Threads.unlock(l)
        # println(" ", i)
    end
end

# Define an example space and optimize it
A = [[float(i), float(j), float(k)] for i = 0:255, j = 0:255, k = 0:255]
@time optimizeSpace!(A)
A[60, 70, 23]

# Run the benchmarking
# Define model
rgb2cam = CAM16Model()

# Define an example point
using Optim
t = [6.0, 4.0, 0.0]
dist(t)
Optim.minimizer(optimize(dist, copy(t)))
Optim.minimizer(optimize(dist, copy(t .+ 0.1), Newton(); autodiff = :forward))
Optim.minimizer(optimize(dist, copy(t .+ 0.1), LBFGS(); autodiff = :forward ))

@benchmark Optim.minimizer(optimize(dist, $copy(t)))

# Define the distance and example configurations that can be reused
dist = distance(t[1], t[2], t[3], 2.0, rgb2cam)
result = DiffResults.HessianResult(copy(t))
config = ForwardDiff.HessianConfig(dist, result, copy(t), ForwardDiff.Chunk{3}())
dist([1.0, 1.0, 0.0])
@benchmark optimize!($copy(t), $rgb2cam, $config, $result)
@profiler for i = 1:100000 optimize!(copy(t), rgb2cam, config, result) end
optimize!(t, rgb2cam, config, result)
dist([2.9999999999999996, 2.9999999999999996, 0.0])

dist([0.38125, 0.38125, 0.006249999999999997])
dist([3.8702555326162194, 3.854765410669157, 0.0])
dist(t)
dist(t)
dist(t)
dist([6.0, 6.0, 0.0])
## TEST FIELD
1.41 *norm(rgb2cam(99, 99, 99) .- rgb2cam(95, 94, 97))^(0.63 * 0.5)


using Plots
using Interact

function plotme(x_value, y_value, z_value, rg, ep, axis)
    rgb2cam = CAM16Model()

    # Define an example point
    t = [x_value, y_value, z_value]

    # Define the distance and example configurations that can be reused
    dist = distance(99.0, 99.0, 99.0, ep, rgb2cam)

    size = 1000

    top_x = x_value + rg
    bottom_x = x_value - rg
    x = bottom_x:(top_x-bottom_x)/size:top_x

    top_y = y_value + rg
    bottom_y = y_value - rg
    y = bottom_y:(top_y-bottom_y)/size:top_y

    function filtered_dist(i, j)
        if axis == 1
            d = dist((x_value, i, j))
        elseif axis == 2
            d = dist((i, y_value, j))
        elseif axis == 3
            d = dist((i, j, z_value))
        end
        return d == Inf ? NaN : d
    end

    z = []
    m = Inf
    argm = [0, 0]
    for i = 1:size
        for j = 1:size
            f = filtered_dist(x[i], y[j])
            if !isnan(f) && f < m
                m = f
                argm = [x[i], y[j]]
            end
            append!(z, filtered_dist(x[i], y[j]))
        end
    end
    println(argm)
    println(m)
    X = repeat(reshape(x, 1, :), length(y), 1)
    Y = repeat(y, 1, length(x))

    return heatmap(x, y, filtered_dist)
end

x_value=widget(0.0:1.0:255.0, label="x_value")
y_value=widget(0.0:1.0:255.0, label="y_value")
z_value=widget(0.0:1.0:255.0, label="z_value")
rg=widget(1:1:10, label="range")
ep=widget(1.0:0.5:10.0, label="eps")
axis = widget(Dict("x" => 1, "y" => 2, "z" => 3), label="axis");

interactive_plot = map(plotme, x_value, y_value, z_value, rg, ep, axis)
plotme(4.0, 4.0, 0.0, 1.0, 2.0, 3)
using Mux, WebIO

dist([3.634, 4.094, 0.0])
dist([4.095959335858323, 4.095959335858323, 0.0])
function app(req) # req is a Mux request dictionary
    vbox(
        hbox(x_value, y_value, z_value, rg, ep, axis), # stack horizontally
        interactive_plot)
end
webio_serve(page("/", app), 8001)
plotly()
gr()

dist = distance(254., 230., 254., 2.0, rgb2cam)
sizez = 1000
top = 222
bottom = 222
x = bottom:(top-bottom)/sizez:top
y = bottom:(top-bottom)/sizez:top

function hello(i, j)
    d = dist((i, j, 99))
    return d == Inf ? NaN : d
end

z = []
for i = 1:sizez
    for j = 1:sizez
        g = hello(x[i], y[j])
        if g != 10000
            append!(z, g)
        end
    end
end
min(z...)
max(z...)
mean(z)

X = repeat(reshape(x, 1, :), length(y), 1)
Y = repeat(y, 1, length(x))

plt = heatmap(x, y, hello)
