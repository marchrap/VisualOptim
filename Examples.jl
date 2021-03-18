using Pkg

Pkg.activate(".")

using VisualOptim
using Images
using JLD

B = [RGB(i/255, j/255, k/255) for i = 0:50, j = 0:50, k = 0:50]
model = CAM16Model()
@profiler optimize_space!(B, 3.0, model)
save("spaces/3.jld", "B", B)
