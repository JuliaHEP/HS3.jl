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

labels, root_vals, errors = read_values("/net/e4-nfs-home.e4.physik.tu-dortmund.de/home/rpelkner/Documents/higgs_analyses/zz/zz_Best_Fit.txt")

using JSON3
nt = JSON3.read("higgs_noES_zzCI.json")
nt = NamedTuple(nt)
labels[1] = "ATLAS_EFF_SIG_4l"
length(labels)
#delete(nt, )
println(Symbol.(labels))
values_julia3 = []
starting_values = []
minus_errors = []
plus_errors = []
for (k, v)  in zip(keys(mode1.result), mode1.result)
    if String(k) in labels
        push!(values_julia3, v)
        push!(minus_errors, (nt)[k][1])
        push!(plus_errors, (nt)[k][2])
    end
end
minus_errors
plus_errors
values_julia3
#for (k, v)  in zip(keys(all_points), all_points)
#    if String(k) in labels
#        push!(starting_values, v)
#    end
#end
#minus_errors = [val[1] for val in values(nt)]
#plus_errors = [val[2] for val in values(nt)]
#values_julia3
# Plot
using Plots
p = plot(size= (1000, 1000))
p = scatter!(root_vals, 1:length(root_vals),  xerr=errors, label=false, ma=4, msw= 0,ms=0, lw =10,coulor=:lightgrey, linecolor=:lightgrey, markercolor=:lightgrey)
p = scatter!(root_vals, 1:length(root_vals), ma=4, msw= 0,ms=12, lw =10,coulor=:lightgrey, linecolor=:lightgrey, markercolor=:lightgrey, label="ROOT Fit Vals & Uncer.", legend=true,)
p = scatter!(values_julia3, 1:length(root_vals), xerr = (values_julia3 - minus_errors, plus_errors - values_julia3), lw =2, label="Julia Fit Vals (LBFGS) & Uncer.", legend=true, markershape=:xcross, colour =:dodgerblue, markercolor=:dodgerblue, linecolor=:dodgerblue, alpha = 0.78)
(800, 600)

#p = scatter!(values_julia3, collect(1:length(root_vals)) .-0.2, xerr = (values_julia3 - minus_errors, plus_errors - values_julia3), label="Julia Fit Vals, LBFGS", legend=true,  markershape=:xcross, ms=6, colour =:green, markercolor=:green, linecolor=:green, alpha = 0.78)

#scatter!(starting_values, 1:length(values), label="Julia Starting Values", legend=true, style="--")
p =plot!(legend=:outerbottom, legendfontsize=10)
p =yticks!(1:length(root_vals), labels, tickfontsize=10)
p =savefig("zz_very_best_fit6.pdf")
collect(1:length(root_vals)) .-0.2



l = Dict(:a => ("k", 2, 3), :b => ("m", 7, 3))
x, y, z, = l

a = [1, 2, 3]
Tuple(a)
using Distributions



#differences = values .- values_julia
#
#p1 = scatter(values, labels, xerr=errors, label=false, legend=false,  title="Values from data1.txt")
#p2 = scatter(differences, 1:length(differences),  label=false, legend=false, title="Difference between data sets")
#
#combined_plot = plot(p1, p2, layout=(1,2), link=:y)
#yticks!(1:length(values), labels)