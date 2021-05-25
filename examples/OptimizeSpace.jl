using Pkg

Pkg.activate(".")

using Revise
using VisualOptim
using Images
using JLD
using Colors

# Get the array of epsilons
epsilons = Array(1.0:5.0)

# Optimize each space using the specified model and save the space
for eps in epsilons
    B = [RGB(i/255, j/255, k/255) for i = 0:255, j = 0:255, k = 0:255]
    model = CAM16Model()
    optimize_space!(B, eps, model, :yes)
    save("spaces/$eps.jld", "B", B)
end
