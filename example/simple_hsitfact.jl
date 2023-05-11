using HS3
using BenchmarkTools
using ValueShapes
example_dict = file_to_dict("./example/sample_normsys.json")

functional_specs  = HS3.make_distspecs(example_dict[:distributions])
data_specs = HS3.make_dataspecs(example_dict[:data])
ll_spec = HS3.generate_likelihoodspec(example_dict[:likelihoods][1])
likelihood = HS3.make_likelihood(ll_spec, functional_specs, data_specs)

using DensityInterface
logdensityof(likelihood, (mu = 200., theta=7, ))
using BAT
using LiteHF, Distributions
a = HS3.make_histfact(functional_specs[2], :mychannel)
a.channel[1]
shape = NamedTupleShape(mu = ScalarShape{Real}(), theta = ScalarShape{Real}())
prior = [LiteHF.FlatPrior(0, 5), Normal()]
prior_dist = product_distribution(prior)
p= shape(prior_dist)
prior = NamedTupleDist(
    a = [Weibull(1.1, 5000), Weibull(1.1, 5000)],
    sigma = Weibull(1.2, 2)
)
unshaped(prior)
(prior)
posterior = PosteriorMeasure(likelihood, shape(prior_dist))
s = bat_sample(posterior, MCMCSampling(mcalg = MetropolisHastings(), nsteps = 10^5, nchains = 4)).result
bat_report(s)
bat_findmode(posterior).result
l = LiteHF.FlatPrior(0, 5)
dump(l)
Distributions.quantile(d::LiteHF.FlatPrior, p::Float64) = d.a + p * (d.b - d.a)
Distributions.params(d::LiteHF.FlatPrior) = (d.a, d.b)



pydict = LiteHF.load_pyhfjson("./example/sample_normsys_lite.json")
pyhfmodel = build_pyhf(pydict)
LL = pyhf_logjointof(pyhfmodel)
LL([200., 7])
priors(pyhfmodel)
ll =  DensityInterface.logfuncdensity(
    function(params)
        return LL(([convert(Float64, params[x]) for x in keys(priors(pyhfmodel))]))
    end)

f((mu =2, theta= 3))
logdensityof(ll, (mu = 1., theta=.1, ))
DensityKind(ll)
post = PosteriorMeasure(ll, shape(prior_dist))
bat_findmode(post).result
###########################################################################################################










a
a  = [LiteHF.FlatPrior(0, 5), Normal()]
a
a = product_distribution(a)
a
#function Distributions.ProductDistribution(dists::AbstractArray{<:Distribution{ArrayLikeVariate{M}},N}) where {M,N}
#    return ProductDistribution{M + N,M,typeof(dists)}(dists)
#end
k = NamedTupleShape(mu = ScalarShape{Real}(), theta = ScalarShape{Real}())
x = k(a)
valshape(x)
unshaped(x)

valshape(k(a)).αs

r = ReshapedDist(a, k)
r.αs
ValueShapes.default_unshaped_eltype(r)
convert(Vector, prior)
prior = HS3.make_histfact(histfact_spec, :mychannel).prior
NamedTupleShape()
using ValueShapes
prior = NamedTupleDist(αs = prior, b = Weibull(1, 2))
varshape(αs)
a

prior = NamedTupleDist(αs = a,)
DensityKind(l)
posterior = PosteriorMeasure(l, k(a))
using Distributions, LiteHF

using Random
LiteHF.rand!(rng::Random.AbstractRNG, d::LiteHF.FlatPrior) = rand(rng, Uniform(d.a, d.b))
s = bat_sample(posterior, MCMCSampling(mcalg = MetropolisHastings(), nsteps = 10^5, nchains = 4)).result
bat_report(s)
@show bat_findmode(posterior).result
l([1, 2])
a.samples[1].modifiers
t = HS3.make_histfact(a, "channel1")
typeof(t) <: HS3.HistfactPDF
p = ( a= 2, b = 3, c= 4)
begin params ->
[params[x] for x in keys(p)]
end

global_unique = Tuple(unique!([:mu, :theta, :mu]))
input_modifiers = [t[2][k] for k in global_unique]
#global_unique = reduce(vcat, [sample[2].modifier_names for (name, C) in channels for sample in C]) |> unique 
p = LiteHF._prior.(input_modifiers)
priors = NamedTuple{global_unique}(p)
prior_log = LiteHF.pyhf_logpriorof(priors)
prior_log([1, 2])
typeof(t[1])
using LiteHF
a = LiteHF.build_channel(t; misc=Dict())
k(a)
v = ("channel1" => a)
name, samples = v
k = LiteHF.build_pyhfchannel(v, [:mu, :theta, :mu])
m = LiteHF.pyhf_loglikelihoodof(k[1], [1, 2.0, 3, 9.0])
m([1., 454.])

s = a.samples[1]
d = Dict(key=>getfield(s, key) for key ∈ fieldnames(HS3.HistFactorySampleSpec))
d[:modifiers] = [m] 
d[:name] = "sample1"
d
chan = Dict()
chan[:samples] = [d]
using LiteHF
l = LiteHF.build_channel(chan; misc=Dict())
l["sample1"].modifier_names
LiteHF._prior(l["sample1"].modifiers[1])
chan
typeof(d[:modifiers])
m = Dict(key=>getfield(m, key) for key ∈ fieldnames(typeof(m)))
m = merge(m, Dict(:name => "mu", :type => "normfactor", :data => nothing))  

HS3.make_histfact(a)
[Dict(:name => "a", :modifiers => [Dict(:name => "mod1", :type => "normfactor") ]) ]
l = Dict(:samples => [Dict(:data => [1, 2, 5], :name => "a", :modifiers => [Dict(:name => "mod1", :type => "shapesys", :data => [2, 3, 5]) ]) ] )
using LiteHF
misc = Dict(:measurements => "mes", :observations => "obs")
k =LiteHF.build_channel(l; misc)

k["a"]
LiteHF.(k)([1, 2, 3])
exp = LiteHF._expkernel(k["a"].modifiers, k["a"].nominal, [1, 2, 3])
LiteHF.pyhf_loglikelihoodof(exp, [1., 2,3])(2)

l[:samples][1][:data]
g = Dict("sample1" => (k, k))
k, v = g
k
LiteHF.build_pyhfchannel(g, Vector())