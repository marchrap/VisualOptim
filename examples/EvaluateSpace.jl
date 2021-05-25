using Pkg

Pkg.activate(".")

using Revise
using VisualOptim
using JLD
using StatsBase
using Images
using Statistics
using Plots
using LinearAlgebra

# Load spaces
optimized = JLD.load("spaces/3.0.jld", "B")
optimized2 = JLD.load("spaces/ipopt/3.0-4-7.jld", "space")
VisualOptim.power(optimized[100, 100, 2].r, optimized[100, 100, 2].g, optimized[100, 100, 2].b)
VisualOptim.power(optimized2[100, 100, 5].r*255, optimized2[100, 100, 5].g*255, optimized2[100, 100, 5].b*255)

# Check power comparisons between two spaces
original_power = 0
new_power = 0
initial_power = 0
for i = 1:256, j = 1:256, k = 1:4
    point = optimized[i, j, k]
    new_power += VisualOptim.power(point.r, point.g, point.b)

    point = optimized2[i, j, k + 3]
    original_power += VisualOptim.power(point.r*255, point.g*255, point.b*255)

    initial_power += VisualOptim.power(i-1, j-1, k + 2)
end
(new_power - original_power)/initial_power

# Evaluate the space for power and color changes
original, modified = VisualOptim.evaluate_space_power(optimized)
println("Power gains per whole space are equal to $(abs(original-modified)/original)")
color_changes, power_changes = VisualOptim.evaluate_space_RGB_change(optimized)

# Obtain some statistics
println("\nR statistics\n-----------------------")
describe(color_changes[:, 1])
println("Standard devation: $(std(color_changes[:, 1]))")

println("\nG statistics\n-----------------------")
describe(color_changes[:, 2])
println("Standard devation: $(std(color_changes[:, 2]))")

println("\nB statistics\n-----------------------")
describe(color_changes[:, 3])
println("Standard devation: $(std(color_changes[:, 3]))")

# Plot the histograms of the space color differences
Plots.histogram(Any[color_changes[:,1], color_changes[:, 2], color_changes[:, 3]])
Plots.histogram(color_changes[:, 3])
Plots.histogram2d(color_changes[:,2], color_changes[:, 3], bins=100)

# Evaluation of the density of the power gains
# Evaluating the different color changes
u = unique(color_changes, dims=1)
d = Dict([(vec, 0) for vec in eachrow(u)])
for i in eachrow(color_changes)
    d[i] += 1
end

# Evaluating how much each discrete change vector contributes to the overall
# power gains
u = unique(color_changes, dims=1)
d = Dict([(vec, 0.) for vec in eachrow(u)])
for (index, i) in enumerate(eachrow(color_changes))
    d[i] += power_changes[index]
end

# Plot the gains
sizes = map(x->x.second, sort(collect(d), by=x->x[2]))
sum(sizes[1:10000])/sum(sizes)
Plots.bar(sizes)

# Obtain the CDF and plot it
distribution = [sizes[1]/sum(sizes)]
for point in sizes[2:end]
    append!(distribution, point/sum(sizes) + distribution[end])
end
Plots.plot(distribution)
Plots.xlabel!("Index")
Plots.ylabel!("Power contribution", fontfamily="Computer Modern", legend=false)
Plots.savefig("../../final_report/figures/compression.pdf")

# This snippet looks at how well the local space can be represented by single
# change vector
kernel_size = 32
local_color_changes = zeros(kernel_size^3, 3)
for i1 = 5:5, j1 = 5:5, k1 = 5:5
    for i2 = 0:kernel_size-1, j2 = 0:kernel_size-1, k2 = 0:kernel_size-1
        i = kernel_size*i1 + i2
        j = kernel_size*j1 + j2
        k = kernel_size*k1 + k2
        index = i*256^2 + j*256 + k + 1
        index2 = i2*kernel_size^2 + j2*kernel_size + k2 + 1
        local_color_changes[index2, 1] = color_changes[index, 1]
        local_color_changes[index2, 2] = color_changes[index, 2]
        local_color_changes[index2, 3] = color_changes[index, 3]
    end
    println("\nR statistics\n-----------------------")
    describe(local_color_changes[:, 1])
    println("Standard devation: $(std(local_color_changes[:, 1]))")

    println("\nG statistics\n-----------------------")
    describe(local_color_changes[:, 2])
    println("Standard devation: $(std(local_color_changes[:, 2]))")

    println("\nB statistics\n-----------------------")
    describe(local_color_changes[:, 3])
    println("Standard devation: $(std(local_color_changes[:, 3]))")
end
Plots.histogram(local_color_changes[:, 1])
savefig("histogram_r_middle_space.pdf")

# Plot the 2D histogram of the changes in a given dimension for different color
# values
local_color_changes = zeros(256*256*256, 2)
d = [Dict() for i = 0:255]
for i = 0:255, j = 0:255, k = 0:255
    index = i*256^2 + j*256 + k+1
    # Remember to change the index here and the name, e.g. 1, i; 2, j; 3, k
    local_color_changes[index, 1] = color_changes[index, 1]
    local_color_changes[index, 2] = i
    dict = d[j+1]
    if haskey(dict, color_changes[index, 2])
        dict[color_changes[index, 2]] += 1
    else
        dict[color_changes[index, 2]] = 1
    end
end

x = []
y = []
z = []
for i = 0:255
    for key in keys(d[i+1])
        append!(x, i)
        append!(y, key)
        append!(z, d[i+1][key])
    end
end
plot(x, y, z)

Plots.histogram2d(local_color_changes[:, 2], local_color_changes[:, 1], show_empty_bins=true)
savefig("histogram_r_whole_space.pdf")

# Read a point from the optimized array
optimized[33, 193, 143]*255

# Check the magnitude of changes on the diagonal and plot them
changes = []
for i = 0:255
    change = norm(i .- [optimized[i+1].r optimized[i+1].g optimized[i+1].b])
    append!(changes, change)
end
Plots.plot(changes)
Plots.xlabel!("Value of the RGB component")
Plots.ylabel!("Magnitude of the change", fontfamily="Computer Modern", legend=false)
Plots.savefig("../../final_report/figures/magnitude_of_change.pdf")
