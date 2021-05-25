using Pkg

Pkg.activate(".")

using Revise
using VisualOptim
using JLD
using StatsBase
using Images
using Statistics
using Plots

# Load space
optimized = JLD.load("spaces/optimized_3.0.jld", "space")

# Load images
img = Images.load("../src/monarch2.png")
Images.save("evaluations/original.png", img)

# If we have integer optimization
img2 = deepcopy(img)
for (index, pixel) in enumerate(img2)
    new = optimized[pixel]
    img2[index] = RGB(new.r/255, new.g/255, new.b/255)
end

# If we have other optimizations
#img2 = VisualOptim.convert_image(img, optimized)

# Save results
Images.save("evaluations/changed.png", img2)

# Evaluate power difference
difference = VisualOptim.compare_power(img, img2)
println("The power gain is $difference")

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
