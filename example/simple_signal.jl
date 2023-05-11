using HS3
using BenchmarkTools


a = Dict(:vals => [0.0, 2.0, 3])
typeof(a) <:AbstractDict{Symbol, <:AbstractArray}
b = Dict( :hi => Dict( :contents => [ -2.5, -3.1 ] ), :lo =>  Dict( :contents => [
    2.2, 3.7 ] ) )
a[:vals]
haskey(a, :vals)
typeof(NamedTuple(a)) <: NamedTuple{(:vals,), <:Tuple{AbstractArray{<:Number}}}

#HS3.generate_all_specs("./example/test.json")

example_dict = file_to_dict("./example/simple_signal.json")
dist_specs = HS3.make_distspecs(example_dict[:distributions])
data_specs = HS3.make_dataspecs(example_dict[:data])
domain_specs = HS3.make_domainspecs(example_dict[:domains])
point_specs = HS3.make_parameterpointspecs(example_dict[:parameter_points])
likelihood_specs = HS3.make_likelihood_specs(example_dict[:likelihoods])
analyses_specs = HS3.make_analysesspecs(example_dict[:analyses])

specs = (analyses = analyses_specs, likelihoods = likelihood_specs, parameter_points = point_specs, domains = domain_specs,
    data = data_specs, functionals = dist_specs,)

ana = HS3.make_analyses(analyses_specs[1], specs)
ll = ana.likelihood
using DensityInterface
@btime logdensityof(ll, (param_mean = 4, param_sigma = 0.5, coef_sig = 0.5, coef_bkg = 0.5,))
prior = ana.prior

data = HS3.make_data(data_specs[1])
using Distributions 
d = MixtureModel([Exponential(1), Normal(params.param_mean, params.param_sigma)], [params.coef_sig, params.coef_bkg])
likelihood = let data = data
    bin_edges = data.edges[1]
    bin_edges_left = bin_edges[1:end-1]
    bin_edges_right = bin_edges[2:end]
    bin_widths = bin_edges_right - bin_edges_left
    bin_centers = (bin_edges_right + bin_edges_left) / 2 
    logfuncdensity(function (params)
    lpdf = 0.0
    idx = eachindex(data.weights)
    for  i in 1:length(data.weights)
         #bin_centers[i] = round(bin_centers[i],  digits=0) 
         expected_counts = bin_widths[i] * Distributions.pdf(MixtureModel([Exponential(1), Normal(params.param_mean, params.param_sigma)], [params.coef_sig, params.coef_bkg]), bin_centers[i])
         #println(bin_centers)
         lpdf += Distributions.logpdf(Distributions.Poisson(expected_counts), data.weights[i])
    end
    #lpdf = Distributions.logpdf(dist(params), data.weights)
     return lpdf
    end)
end
#### usage
@btime logdensityof(likelihood, (param_mean = 4, param_sigma = 0.5, coef_sig = 0.5, coef_bkg = 0.5,))
using BAT
posterior = PosteriorMeasure(ll, prior)
samples = bat_sample(posterior, MCMCSampling(mcalg = MetropolisHastings(), nsteps = 10^6, nchains = 4)).result
using Distributions
bat_report(samples)
using Plots
plot(samples, mean = true, std = true, globalmode = true, marginalmode = true)
println("Mode: $(mode(samples))")
println("Mean: $(mean(samples))")
println("Stddev: $(std(samples))")

a = []
a = push!(a, [])
a[2] = push!(a[2], 5)
a
import Plots
Plots.plot(hist)
Plots.plot!(W)
using Distributions
k = MixtureModel([Normal(4, 0.5), Exponential(1)], [0.5, 0.5])
using StatsBase
hist = append!(Histogram(-0.1:0.1:7), rand(k, 1000000))
print((hist.weights))
samples = bat_sample(k, MCMCSampling(mcalg = MetropolisHastings(), nsteps = 10^5, nchains = 4)).result
