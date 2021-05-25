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

# Define the custom color struct
# Necessary as the RGB one is immutable
mutable struct custom_color
    r::Float64
    g::Float64
    b::Float64
end

# Space optimization function with a progress meter
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

# Load the required images
img = load("../src/fruits.png")
img2 = load("../src/monarch2.png")
println(length(countmap(img)))
println(length(countmap(img2)))

# Obtain the unique RGB triplets that are saved as spaces
total = vec(img)
total = unique!(total)
append!(total, vec(img2))
total = unique!(total)
save("spaces/space.jld", "space", total)

# Load and optimize the spaces - i.e. here the unique pixels of an image
original = load("spaces/optimized_5.0.jld", "space")
space = [custom_color(pixel.r*255, pixel.g*255, pixel.b*255) for pixel in original]
@time optimize_space_int!(space, epss)

# Save the results in a form of a dictionary - this is better when it comes
# to accessing the results that the other possible datatypes
a = Dict()
for (index, value) in enumerate(space)
    a[original[index]] = value
end
save("spaces/optimized_$epss.jld", "space", a)

# For all the epsilons check the power gains for all the images and save
# them so that they can be put into a grid
for eps = 1.:5.
    # Load space
    optimized = JLD.load("spaces/optimized_$eps.jld", "space")

    # Evaluate images
    for i in ["fruits.png", "monarch2.png"]
        # Load images, convert img2 and save results
        img = load("../src/$i")
        img2 = deepcopy(img)
        for (index, pixel) in enumerate(img2)
            new = optimized[pixel]
            img2[index] = RGB(new.r/255, new.g/255, new.b/255)
        end
        diff = VisualOptim.compare_power(img, img2)
        save("../../final_report/figures/images/CIECAM16_IPOPT/$(@sprintf("%.4f",diff))_$(eps)_$i", img2)
    end
end
