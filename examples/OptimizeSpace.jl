using Pkg

Pkg.activate(".")

using VisualOptim
using Images
using JLD

epsilons = [1.0, 2.0, 3.0, 4.0, 5.0]

for eps in epsilons
    B = [RGB(i/255, j/255, k/255) for i = 0:255, j = 0:255, k = 0:255]
    model = CAM16Model()
    optimize_space!(B, eps, model)
    save("spaces/$eps.jld", "B", B)
end
