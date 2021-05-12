# Necessary imports
using LinearAlgebra
using Colors
using Images
using CSV
using DataFrames
using BenchmarkTools
using StaticArrays
using Roots
using Printf

"""
The overall loop that obtaines the optimal converted space given a max
deviation eps.

Params:
    eps: The maximum allowable deviation from the original point.
"""
function optimize(eps)
    # Define the model parameters
    Dα = SA[5.12040476e-02 7.12853167e-03 -9.22189345e-03;
            7.12853167e-03 1.08748289e-03 -2.55714768e-04;
            -9.22189345e-03 -2.55714768e-04 4.61453920e-03]

    p = SA[-1.25929564e+00; -1.25505864e-01; 1.09027499e-01]

    # Eigenvalue decomposition of the matrix
    eigen_values, eigen_vectors = eigen(Dα)

    # Define the required array
    optimized = Array{RGB}(undef, 256, 256, 256)

    # Spread the work between threads
    Δ = 256 / Threads.nthreads()
    indexes = [(Int(Δ*(i-1)), Int(Δ*i-1)) for i in 1:Threads.nthreads()]

    # Run the loop of convert -> optimize -> write to array
    Threads.@threads for index in indexes
        task(optimized, index, eigen_values, eigen_vectors, Dα, p, eps)
    end

    # Return the final array
    return optimized
end

function task(results::Array{RGB, 3}, indexes::Tuple{Int64, Int64}, eigen_values::SArray{Tuple{3},Float64,1,3}, eigen_vectors::SArray{Tuple{3,3},Float64,2,9}, power_matrix::SArray{Tuple{3,3},Float64,2,9}, p::SArray{Tuple{3},Float64,1,3}, eps::Float64)
    #=
        The function that is run by the threads and simply runs a for loop that
        iterates over the given range and solves the problem.

        Params:
            results: The array of results to which the solutions should be placed.
            indexes: A tuple that containes the range of indexes that this routine
                     should consider solving.
            eigen_values: The eigen_values as computed by optimize.
            eigen_vectors: The eigen_vectors as computed by optimize.
            power_matrix: The power matrix as defined in optimize.
            p: The p as defined in optimize (the linear term of power).
            eps: Maximum allowable deviation from the original point.
    =#

    # Spread the indexes
    min = indexes[1]
    max = indexes[2]

    # Run the for loop to evaluate the points
    for r in min:max, g in 0:255, b in 0:255
        solve(results, (r, g, b), eigen_values, eigen_vectors, power_matrix, p, eps)
    end
end

function solve(results::Array{RGB, 3}, index::Tuple{Int64, Int64, Int64}, eigen_values::SArray{Tuple{3},Float64,1,3}, eigen_vectors::SArray{Tuple{3,3},Float64,2,9}, power_matrix::SArray{Tuple{3,3},Float64,2,9}, p::SArray{Tuple{3},Float64,1,3}, eps::Float64)
    #=
        The function that is solving the problem using a simple line search.

        Params:
            results: The array of results to which the solutions should be placed.
            indexes: A tuple that containes the range of indexes that this routine
                     should consider solving.
            eigen_values: The eigen_values as computed by optimize.
            eigen_vectors: The eigen_vectors as computed by optimize.
            power_matrix: The power matrix as defined in optimize.
            p: The p as defined in optimize (the linear term of power).
            eps: Maximum allowable deviation from the original point.
    =#

    # Convert index to colorspace
    v = colorspace_forward(index[1], index[2], index[3])

    # Define bs and βs
    b = -1 .* (v' * power_matrix .+ 1/2 .* p')

    # Define β
    β = eigen_vectors' * b'

    # Obtain the optimal value of vectors
    if β[1] ≈ 0
        # Degenerate case
        temp = eigen_values[2:end] .- eigen_vectors[1]
        c = β[2:end] ./ temp
        u = eigen_vectors[:, 2:end] * c .+ v
    else
        # Non-degenerate case
        evaluator = function_generator(β, eigen_values, eps)
        μ_l = abs(β[1])/eps-eigen_values[1]
        μ_u = norm(b)/eps-eigen_values[1]
        μ = (μ_l+μ_u)/2
        step = 1.0

        # Newton's method
        while abs(step) > 1e-8
            step = evaluator(μ)
            μ = μ - step
        end

        # Obtain optimal u
        temp = eigen_values .+ μ
        c = β ./ temp
        u = eigen_vectors * c .+ v
    end

    # Convert from back to RGB and change the results matrix
    converted = colorspace_backward(u[1], u[2], u[3])
    @inbounds results[index[1]+1, index[2]+1, index[3]+1] = converted
end

function function_generator(β, λ, eps)
    #=
        A function that generates a function that gives the change required by
        the Halley's method.

        Params:
            β: Beta as defined in the problem statement.
            λ: Lambda as defined in the problem statement.
            eps: The maximum allowable deviation from the optimal point.
    =#

    a = function func(μ::Float64)
        #=
            Implements the Halley's method and its delta x

            Params:
                μ: The point at which current method is evaluated
        =#
        temp = (λ .+ μ)
        c = (β.^2) ./ (temp.^2)
        c_2 = c ./ temp
        c_3 = c_2 ./ temp
        f = sum(c)-eps^2
        fp = -2 * sum(c_2)
        fpp = 6* sum(c_3)
        return f/fp/(1-f/fp*fpp/2/fp)
    end

    return a
end

function colorspace_forward(r::Int64, g::Int64, b::Int64)
    #=
        Converts the RGB colors to the different color spaces. Allows for different
        colorspaces.

        Params:
            r: The red part of the color.
            g: The green part of the color.
            b: The blue part of the color.
    =#
    color = convert(Lab, RGB(r/255, g/255, b/255))
    return SA[color.l, color.a, color.b]
end

function colorspace_backward(x::Float64, y::Float64, z::Float64)
    #=
        Converts the colors to the RGB color space. Allows for different
        colorspaces.

        Params:
            x: The first part of the color.
            y: The second part of the color.
            z: The thirtd part of the color.
    =#
    # Allows for different colorspaces
    color = convert(RGB, Lab(x, y, z))
    return color
end

function convert_pixel(pixel::RGB, optimized::Array{RGB, 3})
    #=
        A helper function that allows to read the the result from the optimized
        array and protects against the NaN values.

        Params:
            pixel: The RGB triplet that specifies which color is being optimized.
            optimized: The optimized array.
    =#
    # Obtain the value
    value = optimized[convert(Int64, 255*pixel.r)+1, convert(Int64, 255*pixel.g)+1, convert(Int64, 255*pixel.b)+1]

    # Check if the value is NaN
    if isnan(value.r)
        return pixel
    end
    return value
end

function evaluate_power(color::RGB, df::DataFrame)
    #=
        A helper function that allows to evaluate the power of a given color.

        Params:
            color: An RGB triplet that specifies the color.
            df: The dataframe corresponding to the measured power data.
    =#
    # Convert the colors to ints
    R = convert(Int64, round(color.r*255))
    G = convert(Int64, round(color.g*255))
    B = convert(Int64, round(color.b*255))

    # Read the dataframe
    total = 0
    if R != 0
        total += df[(df.R .== R) .& (df.G .== 0) .& (df.B .== 0), :currentNow][1]
    end
    if G != 0
        total += df[(df.R .== 0) .& (df.G .== G) .& (df.B .== 0), :currentNow][1]
    end
    if B != 0
        total += df[(df.R .== 0) .& (df.G .== 0) .& (df.B .== B), :currentNow][1]
    end

    return total
end

function compare_power(img1, img2)
    #=
        Compares the power of two images.

        Params:
            img1: The first image as loaded by Julia.
            img2: The second image as loaded by Julia.
    =#
    original_power = 0
    new_power = 0

    # Read the power data
    power_data = CSV.read("project/optimisation/power_data/helsinki1.csv")
    power_data= select(power_data, [:currentNow, :R, :G, :B])
    power_data[:, :currentNow] = - power_data[:, :currentNow]

    # Evaluate the power
    for i in 1:size(img1)[1], j in 1:size(img1)[2]
        original_power += evaluate_power(img1[i, j], power_data)
        new_power += evaluate_power(img2[i, j], power_data)
    end

    # Print the results
    println("Original power $original_power; New power $new_power")

    return (new_power-original_power)/original_power
end

# Tests
# Benchamrk timing
@benchmark optimize(5.0)
@benchmark optimize(5.0)

# Run the algorithm for some example images
eps = 8.0
optimized = optimize(eps)

for i in ["photo.png", "cat.png", "fruits.png", "monarch.png"]
    img = load("project/optimisation/src/$i")
    b = map(x -> convert_pixel(x, optimized), img)
    diff = compare_power(img, b)
    save("project/optimisation/src/optimized_$(@sprintf("%.4f",diff))_$(eps)_$i", b)
end
