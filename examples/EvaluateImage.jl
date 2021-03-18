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
optimized = JLD.load("spaces/3.jld", "B")

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
