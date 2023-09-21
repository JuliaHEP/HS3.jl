julia_fit = mode1.result



function read_values(filename)
    labels = String[]
    values = Float64[]
    errors = Float64[]
    
    open(filename, "r") do file
        for line in eachline(file)
            m = match(r"(\S+)\s+([-+]?\d*\.\d+e[+-]\d+)\s+\+/-\s+([-+]?\d*\.\d+e[+-]\d+)", line)
            if m !== nothing
                push!(labels, m.captures[1])
                push!(values, parse(Float64, m.captures[2]))
                push!(errors, parse(Float64, m.captures[3]))
            end
        end
    end
    return labels, values, errors
end

labels, root_vals, errors = read_values("zz_Best_Fit.txt")
labels[1] = "ATLAS_EFF_SIG_4l"
labels
values_julia3 = []
starting_values = []
for (k, v)  in zip(keys(mode1.result), mode1.result)
    if String(k) in labels
        push!(values_julia3, v)
    end
end

for (k, v)  in zip(keys(all_points), all_points)
    if String(k) in labels
        push!(starting_values, v)
    end
end
vals == values_julia3
values
values_julia3
all_points
starting_values

# Plot
using Plots
scatter(root_vals, 1:length(root_vals),  xerr=errors, label="ROOT Fit Vals & Uncer.", legend=true, ma=4, coulor=:blue, linecolor=:blue)
scatter!(values_julia3, 1:length(root_vals), label="Julia Fit Vals, LBFGS", legend=true,  markershape=:xcross, ms=6)
#scatter!(starting_values, 1:length(values), label="Julia Starting Values", legend=true, style="--")
plot!(legend=:outerbottom)
yticks!(1:length(root_vals), labels)
savefig("zz_very_best_fit4.pdf")


#differences = values .- values_julia
#
#p1 = scatter(values, labels, xerr=errors, label=false, legend=false,  title="Values from data1.txt")
#p2 = scatter(differences, 1:length(differences),  label=false, legend=false, title="Difference between data sets")
#
#combined_plot = plot(p1, p2, layout=(1,2), link=:y)
#yticks!(1:length(values), labels)