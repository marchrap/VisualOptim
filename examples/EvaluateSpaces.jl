# Allows for checking the distribution of the power gains and changes introduced
# over different epsilons

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

# Define the required parameters
overall_gains = []
overall_color_changes = Int[]
color_group = String[]
eps_group = Int[]

# Iterate over all epsilons
for eps = 1.:5.
    # Print epsilon for debugging purposes
    println(eps)

    # Load spaces
    optimized = JLD.load("spaces/eigenvalue/$eps.jld", "B")

    # Evaluate the space for power and color changes
    # Append the color changes to the different arrays for different colors
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

# Evaluate original total power
total_power = 0
for i = 0:255, j = 0:255, k = 0:255
    total_power += VisualOptim.power(i, j, k)
end
total_power = total_power/256^3

# Gains plot
Plots.plot(overall_gains, fontfamily="Computer Modern", label="barrier")
Plots.scatter!(overall_gains, label=nothing)
xlabel!("ϵ")
ylabel!("power gain", fontfamily="Computer Modern")
Plots.savefig("../../final_report/figures/CIECAM_barrier_comparison.pdf")
# Defined using the data print at each evaluation
a = (total_power .- [155.78, 148.88, 139.65, 129.59, 119.34]) / total_power
Plots.plot!(a, ribbon=0.002794, linecolor=:red, label="IPOPT")
Plots.scatter!(a, label=nothing)

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
    optimized = JLD.load("spaces/optimized_$eps.jld", "B")
    for i in ["photo.png", "cat2.png", "fruits.png", "monarch2.png"]
        img = load("../src/$i")
        b = map(x -> VisualOptim.convert_pixel(x, optimized), img)
        diff = VisualOptim.compare_power(img, b)
        save("../../final_report/figures/images/CIECAM16_IPOPT/$(@sprintf("%.4f",diff))_$(eps)_$i", b)
    end
end
