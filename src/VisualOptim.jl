module VisualOptim

using ProgressMeter
using Optim
using Images
using Plots

include("ColorSpaces.jl")
include("Models.jl")
include("Helpers.jl")

export optimize_space!, CAM16Model, evaluate_space, convert_pixel, histogram

end # module
