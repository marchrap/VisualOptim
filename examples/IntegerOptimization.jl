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

# Run the optimization for the required
init = parse(Int64, ARGS[1])
final = parse(Int64, ARGS[2])
epss = 3.0
space = [custom_color(i, j, k) for i = init:final, j = init:final, k = init:final]
@time optimize_space_int!(space, epss)
save("spaces/ipopt/$epss-$init-$final.jld", "space", space)
