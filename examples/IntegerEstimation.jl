using Pkg

Pkg.activate(".")

using Revise
using VisualOptim
using JuMP
using Juniper
using BenchmarkTools
using Ipopt
using Images
using ProgressMeter
using JLD
using StatsBase

mutable struct custom_color
    r::Float64
    g::Float64
    b::Float64
end

function optimize_space_local!(space, epss, progress_bar, integer=nothing)
    # Initialize the threads with the corresponding model
    # Create the model
    optimizer = Juniper.Optimizer
    nl_solver = optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)
    model = Model(optimizer_with_attributes(optimizer, "nl_solver"=>nl_solver, "log_levels"=>[]))
    color_model = CAM16Model()

    # Initialize the parameters of the model
    @NLparameter(model, x_r == 0)
    @NLparameter(model, x_g == 0)
    @NLparameter(model, x_b == 0)

    # Initialize the distance function
    function du(r, g, b)
        # Obtain the point position of the color in the required space
        x_n, y_n, z_n = color_model(r, g, b)

        # Calculate ΔE
        return 1.41 * ((value(x_r) - x_n)^2 + (value(x_g) - y_n)^2 +
         (value(x_b) - z_n)^2)^(0.63 * 0.5)
    end

    # Register the necessary nonlinear functions
    register(model, :distance, 3, du, autodiff=true)
    register(model, :power, 3, VisualOptim.power, autodiff=true)

    # Create the start variables
    start_point = space[1]
    start_r = convert(Int64, start_point.r)
    start_g = convert(Int64, start_point.g)
    start_b = convert(Int64, start_point.b)

    # Initialize the variables
    if integer == nothing
        @variable(model, 0 <= r <= 255)
        @variable(model, 0 <= g <= 255)
        @variable(model, 0 <= b <= 255)
    else
        @variable(model, 0 <= r <= 255, Int)
        @variable(model, 0 <= g <= 255, Int)
        @variable(model, 0 <= b <= 255, Int)
    end

    # Set the constraints and objectives
    @NLconstraint(model, constraint, distance(r, g, b) <= epss)
    @NLobjective(model, Min, power(r, g, b))
    failed = 0

    # Optimize the model and save the results
    @showprogress for point in space
        # Convert the point to RGB values and then to the used color model
        r_original = convert(Int64, point.r)
        g_original = convert(Int64, point.g)
        b_original = convert(Int64, point.b)
        converted_point = color_model(r_original, g_original, b_original)

        # Set the parameters
        set_value(x_r, converted_point[1])
        set_value(x_g, converted_point[2])
        set_value(x_b, converted_point[3])

        # Set start values
        set_start_value(r, start_r)
        set_start_value(g, start_g)
        set_start_value(b, start_b)

        # Optimize
        try
            optimize!(model)
            point.r = value(r)
            point.g = value(g)
            point.b = value(b)

            # Update the start values
            start_r = value(r)
            start_g = value(g)
            start_b = value(b)
        catch e
            failed += 1

            point.r = r_original
            point.g = g_original
            point.b = b_original

            # Update the start values
            start_r = r_original
            start_g = g_original
            start_b = b_original
        end
    end
    println("$failed failed optimizations")
end

function optimize_space_int!(space, epss, integer=nothing)
    # Create the progress bar
    progress_bar = Progress(length(space), 1, "Computing the matrix... ", 50)
    update!(progress_bar, 0)

    # Run over the points
    optimize_space_local!(space, epss, progress_bar, integer)
end

# Define parameters
epss = parse(Float64, ARGS[1])

# Sample N random points
N = 100000
indexes_R = sample(0:255, N)
indexes_G = sample(0:255, N)
indexes_B = sample(0:255, N)

# Create the required spaces
space = [custom_color(indexes_R[i], indexes_G[i], indexes_B[i]) for i = 1:N]
original_space = deepcopy(space)

# Optimize the space
optimize_space_int!(space, epss)

# Analysis of power
# Define power containers
original_power = 0
new_power = 0

# Define the points
x = Float64[]
y = Float64[]

x_original = Float64[]
y_original = Float64[]
for i = 1:N
    original_power += VisualOptim.power(original_space[i].r, original_space[i].g, original_space[i].b)
    new_power += VisualOptim.power(space[i].r, space[i].g, space[i].b)
    if i % 10 == 0
        append!(x, i)
        append!(y, new_power/i)
    end

    if i % 10 == 0
        append!(x_original, i)
        append!(y_original, original_power/i)
    end
end
println(mean(new_power))
Plots.scatter(x, y, markersize=2)
Plots.scatter!(x_original, y_original, markersize=2)

# Evaluate the standard deviation test for 100 runs each 10^6 points
B = 100
Ni = 1000000
x_original = [0. for i=1:Ni/10]
y_original = [[0. for i=1:Ni/10] for k=1:B]
for j = 1:B
    indexes_Ri = sample(0:255, Ni)
    indexes_Gi = sample(0:255, Ni)
    indexes_Bi = sample(0:255, Ni)
    original_spacei = [custom_color(indexes_Ri[i], indexes_Gi[i], indexes_Bi[i]) for i = 1:Ni]
    original_power = 0

    for i = 1:Ni
        original_power += VisualOptim.power(original_spacei[i].r, original_spacei[i].g, original_spacei[i].b)
        if (i-1) % 10 == 0
            x_original[Int(floor(i/10) + 1)] = i
            y_original[j][Int(floor(i/10) + 1)] = original_power/i
        end
    end
end

# Plot the mean with the standard deviation ribbon
Plots.scatter(x_original, mean(y_original)[1:1000:end], ribbon=std(y_original)[1:1000:end], markersize=2)
Plots.plot(x_original[1:10000], y_original[1][1:10000])

# Plot the standard deviation against the number of iterations
Plots.plot(x_original[1:10000], std(y_original)[1:10000])
Plots.xlabel!("Number of iterations")
Plots.ylabel!("Standard deviation of the estimate", legend=false, fontfamily="Computer Modern")
Plots.savefig("../../final_report/figures/monte_carlo_original.pdf")

# Plot the normalized standard deviation
Plots.plot(x_original[1000:end], (2*std(y_original)[1000:end]/(total_power/256^3)), xaxis=:log)#./mean(y_original))[1000:end])
Plots.xlabel!("Number of iterations")
Plots.ylabel!("2 times the normalized standard deviation", legend=false, fontfamily="Computer Modern")
Plots.savefig("../../final_report/figures/monte_carlo_zoom.pdf")

# Various power estimate accesses 
y_original
mean(y_original)
(2*std(y_original)./mean(y_original))[10000]
x_original[1000]
total_power = 0
for i = 0:255, j = 0:255, k = 0:255
    total_power += VisualOptim.power(i, j, k)
end
original_power/Ni
total_power/256^3
new_power/N - 0.1 * (original_power/N - total_power/256^3)
new_power/N
(new_power - original_power)/original_power
