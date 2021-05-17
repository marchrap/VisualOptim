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
using Printf

overall_gains = []
overall_color_changes = Int[]
color_group = String[]
eps_group = Int[]

for eps = 1.:5.
    println(eps)
    # Load spaces
    optimized = JLD.load("spaces/CIELAB/$eps.jld", "B")

    # Evaluate the space for power and color changes
    original, modified = VisualOptim.evaluate_space_power(optimized)
    println("Power gains per whole space are equal to $(abs(original-modified)/original)")
    append!(overall_gains, abs(original-modified)/original)
    color_changes, power_changes = VisualOptim.evaluate_space_RGB_change(optimized)
    append!(overall_color_changes, color_changes[:, 1])
    append!(color_group, repeat(["R"], inner=length(color_changes[:, 1])))
    append!(overall_color_changes, color_changes[:, 2])
    append!(color_group, repeat(["G"], inner=length(color_changes[:, 2])))
    append!(overall_color_changes, color_changes[:, 3])
    append!(color_group, repeat(["B"], inner=length(color_changes[:, 3])))
    append!(eps_group, repeat([eps], inner=length(color_changes)))
end

# Gains plot
Plots.plot!(overall_gains, fontfamily="Computer Modern", label="barrier")
Plots.scatter!(overall_gains, label=nothing)
xlabel!("ϵ")
ylabel!("power gain", fontfamily="Computer Modern", legend = false)
Plots.savefig("../../final_report/figures/CIELAB_barrier_comparison.pdf")

# Boxplot for different ColorSpaces
indexes = sample(1:length(eps_group), 1000000, replace=false)
groupedviolin(eps_group[indexes], overall_color_changes[indexes], group=color_group[indexes], color=[:royalBlue :green :red], range=1.5)
groupedboxplot!(eps_group[indexes], overall_color_changes[indexes], group=color_group[indexes], fillalpha=0.75, fillcolor=[:lightblue :lightgreen :orange], linecolor=[:black], range=1.5, outliers=false, label=nothing, showmeans=true)
xlabel!("ϵ")
ylabel!("changes for different color values", fontfamily="Computer Modern")
savefig("../../final_report/figures/CIELAB_barrier_change_colors.pdf")

# Evaluate some images
for eps = 1.:5.
    # Load space
    optimized = JLD.load("spaces/CIECAM/$eps.jld", "B")
    for i in ["photo.png", "cat.png", "fruits.png", "monarch.png"]
        img = load("../src/$i")
        b = map(x -> convert_pixel(x, optimized), img)
        diff = VisualOptim.compare_power(img, b)
        save("../../final_report/figures/images/CIECAM16/$(@sprintf("%.4f",diff))_$(eps)_$i", b)
    end
end


for i in ["photo.png", "cat.png", "fruits.png", "monarch.png"]
    img = load("../src/$i")
    println(length(countmap(img)))
end
