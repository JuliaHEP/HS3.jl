using HS3
using DensityInterface
using BenchmarkTools
dict = file_to_dict("./example/test_interp0d.json")
specs = HS3.generate_specs(dict)

func = HS3.make_functional(specs.distributions.model_channel1, HS3.topological_sort(merge(specs.distributions, specs.functions)), 0)
@benchmark func((mu=1, a= 2))
analysis = HS3.make_analyses(specs.analyses[1], specs)
a = logdensityof(analysis.likelihood, (mu =0.2, a=2,))
HS3.topological_sort(merge(specs.distributions, specs.functions))
using DensityInterface
a = logdensityof(analysis.likelihood)
specs.functions.mu_interp

a = (b =1, c=2)
keys(a)
a[:b]

using StatsBase
function generate_likelihood(dist, data::StatsBase.Histogram)  
    bin_edges = data.edges[1]
    bin_edges_left = bin_edges[1:end-1]
    bin_edges_right = bin_edges[2:end]
    bin_widths = bin_edges_right - bin_edges_left
    bin_centers = (bin_edges_right + bin_edges_left) / 2 
    likelihood = function (params)
        lpdf = 0.0
        for i in 1:length(data.weights)
            lpdf += Distributions.logpdf(Distributions.Poisson(bin_widths[i] * Distributions.pdf(dist(params), bin_centers[i])), data.weights[i])
        end
        return lpdf
    end
    return params -> likelihood(params)
end
 
10000000 * 2e-7