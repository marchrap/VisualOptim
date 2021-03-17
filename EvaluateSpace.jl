using Pkg

Pkg.activate(".")

using Revise
using VisualOptim
using JLD
using StatsBase
using Images
using Statistics
using Plots

using StatsPlots
using Clustering

using LinearAlgebra



# Load space
optimized = JLD.load("spaces/3.jld", "B")

# Evaluate the space for power and color changes
original, modified = VisualOptim.evaluate_space_power(optimized)
println("Power gains per whole space are equal to $(abs(original-modified)/original)")
changes = VisualOptim.evaluate_space_RGB_change(optimized)

# Obtain some statistics
println("\nR statistics\n-----------------------")
describe(changes[:, 1])
println("Standard devation: $(std(changes[:, 1]))")

println("\nG statistics\n-----------------------")
describe(changes[:, 2])
println("Standard devation: $(std(changes[:, 2]))")

println("\nB statistics\n-----------------------")
describe(changes[:, 3])
println("Standard devation: $(std(changes[:, 3]))")

# Load images
img = Images.load("../src/monarch.png")
Images.save("evaluations/original.png", img)
img2 = VisualOptim.convert_image(img, optimized)
Images.save("evaluations/changed.png", img2)

# Evaluate power difference
diff = VisualOptim.compare_power(img, img2)
println("The power gain is $diff")

# Evaluate color changes
VisualOptim.compare_colors(img, img2)

# Plot next to each other
VisualOptim.plot_side_by_side(img, img2)
savefig("evaluations/side_by_side.pdf")

# Plot the image difference
VisualOptim.image_compare_histogram(img, img2)
savefig("evaluations/histogram.pdf")

# Plot differences in other colors
VisualOptim.image_compare_histogram_color(img, img2, "R")
savefig("evaluations/histogram_R.pdf")
VisualOptim.image_compare_histogram_color(img, img2, "G")
savefig("evaluations/histogram_G.pdf")
VisualOptim.image_compare_histogram_color(img, img2, "B")
savefig("evaluations/histogram_B.pdf")


for i in ["photo.png", "cat.png", "fruits.png", "monarch.png"]
    img = load("../src/$i")
    R, G, B = VisualOptim.histogram(img)
    plot(plot(img), Plots.histogram(R), Plots.histogram(G), Plots.histogram(B), layout = @layout([ a [b; c; d]]))
    b = map(x -> convert_pixel(x, A), img)
    R, G, B = VisualOptim.histogram(b)
    plot(plot(b), Plots.histogram(R), Plots.histogram(G), Plots.histogram(B), layout = @layout([ a [b; c; d]]))
    diff = compare_power(img, b)
end


sampled = changes[rand(1:length(changes[:, 1]), 100000), :]
Plots.histogram(Any[changes[:,1], changes[:, 2], changes[:, 3]])
Plots.histogram2d(changes[:,2], changes[:, 3], bins=100)
Plots.histogram(changes[:, 3])
cornerplot(changes)
h = fit(Histogram, changes)
h.edges
h.weights

u = unique(changes, dims=1)
d = Dict([(vec, 0) for vec in eachrow(u)])
for i in eachrow(changes)
    d[i] += 1
end

sizes = map(x->x.second, sort(collect(d), by=x->x[2], rev=true))
sum(sizes[1:1000])/sum(sizes)

Plots.bar(sizes)

Plots.bar(map(x->x.second, sort(collect(d), by=x->x[2], rev=true)))
StatsBase.countmap(changes[:, :])
Plots.bar(StatsBase.countmap(changes[:,1]), bar_width=1, alpha=0.3)
Plots.bar!(StatsBase.countmap(changes[:,2]), bar_width=1, alpha=0.3)
Plots.bar!(StatsBase.countmap(changes[:,3]), bar_width=1, alpha=0.3)
optimized[100, 100, 100]
scatter3d(sampled[:, 1], sampled[:, 2], sampled[:, 3], legend=false)
