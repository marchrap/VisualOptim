using Pkg

Pkg.activate(".")

using Revise
using VisualOptim
using Images
using JLD
using Colors

epsilons = [2.0]

for eps in epsilons
    B = [RGB(i/255, j/255, k/255) for i = 0:255, j = 0:255, k = 0:255]
    model = Colors.DE_AB()
    optimize_space!(B, eps, model, :yes)
    save("spaces/$eps.jld", "B", B)
end
