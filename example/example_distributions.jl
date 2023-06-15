using HS3
using BenchmarkTools

@benchmark dict = file_to_dict("./example/a.json")
dict = file_to_dict("./example/a.json")
specs = HS3.generate_specs(dict)
@benchmark HS3.generate_specs(dict)
k = ["a", "b"]
eachindex(k)
dist = specs.distributions.regularization_ttW

@benchmark a = HS3.make_functional(dist, HS3.topological_sort(merge(specs[:distributions], specs[:functions])))
@btime a((mu_prime_2 = 2,
ttW_Bin_001_Ac =1,
ttW_Bin_001_muInc = 1,
ttW_Bin_002_Ac =1,
ttW_Bin_002_muInc = 1,
ttW_Bin_003_Ac = 1,
ttW_Bin_003_muInc = 1,
ttW_Bin_004_Ac = 1,
ttW_Bin_004_muInc = 1,
ttW_Bin_005_Ac = 1,
ttW_Bin_005_muInc = 1,
ttW_Bin_006_Ac = 1,
ttW_Bin_006_muInc = 1,
ttW_Bin_007_Ac = 1,
ttW_Bin_007_muInc = 1,
ttW_Bin_008_Ac = 1,
ttW_Bin_008_muInc = 1,
ttW_Bin_009_Ac = 1,
ttW_Bin_009_muInc = 1,
ttW_Bin_010_Ac = 1.,
ttW_Bin_010_muInc = 1,
ttW_Bin_011_Ac = 1,
ttW_Bin_011_muInc = 1,
ttW_Bin_012_Ac = 1,
ttW_Bin_012_muInc = 1,
ttW_Bin_013_Ac = 1,
ttW_Bin_013_muInc = 1))
a[1]
d = a((mu_prime_2 = 2, mu_prime_1 =4 ))
d[1]
d[2]
using Distributions
k = MvNormal([0., 0], reduce(hcat, [[1, 0.], [0,1]]))
h = Normal()

Distributions.logpdf(h, 1)
#@btime HS3._val_content(Val(:a))
##@btime Symbol("a")
#example_dict = file_to_dict("./example/a.json") ####~ 50ms / 36ms Rastaban
@benchmark example_dict = file_to_dict("./example/a.json") 
#example_dict[:distributions][2].samples[1]
oldstd = stdout
redirect_stdout(Base.DevNull())
Val{:mu_prime_2} in [Val{:mu_prime_2}]
#redirect_stdout(oldstd)
#BenchmarkTools.DEFAULT_PARAMETERS.seconds = 5
#BenchmarkTools.DEFAULT_PARAMETERS.samples = 1000
#BenchmarkTools.DEFAULT_PARAMETERS.evals = 1
#HS3.generate_specs(example_dict) # ~9s
specs = HS3.generate_specs(example_dict) # ~9s, Rastaban: max. 0.5s
@benchmark HS3.generate_specs(example_dict)
fieldnames((a=2, b=3))
specs.likelihoods[2]
@benchmark HS3.make_likelihood(specs.likelihoods[2], merge(specs.distributions, specs.functions), specs.data)
a
using DensityInterface
logdensityof(a, valus)
valus = (
    alpha_ttbb_XS=1,
    ttW_Bin_001_Ac =1,
    ttW_Bin_001_muInc = 1,
    ttW_Bin_002_Ac =1,
    ttW_Bin_002_muInc = 1,
    ttW_Bin_003_Ac = 1,
    ttW_Bin_003_muInc = 1,
    ttW_Bin_004_Ac = 1,
    ttW_Bin_004_muInc = 1,
    ttW_Bin_005_Ac = 1,
    ttW_Bin_005_muInc = 1,
    ttW_Bin_006_Ac = 1,
    ttW_Bin_006_muInc = 1,
    ttW_Bin_007_Ac = 1,
    ttW_Bin_007_muInc = 1,
    ttW_Bin_008_Ac = 1,
    ttW_Bin_008_muInc = 1,
    ttW_Bin_009_Ac = 1,
    ttW_Bin_009_muInc = 1,
    ttW_Bin_010_Ac = 1.,
    ttW_Bin_010_muInc = 1,
    ttW_Bin_011_Ac = 1,
    ttW_Bin_011_muInc = 1,
    ttW_Bin_012_Ac = 1,
    ttW_Bin_012_muInc = 1,
    ttW_Bin_013_Ac = 1,
    ttW_Bin_013_muInc = 1
    )
typeof(valus)
    eh = ((-0.006452301828157481, 0.021704004705953794, -0.013718764463578559, -0.00015046859696732362, -0.004281092937398201, 0.02663982668771725, -0.0042507829075222325, -0.021982735947919663, -0.012532373928357154, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0))
l = [x for x in eh]
typeof(l)
length(Tuple{Val{:mu_prime_2}, Val{:mu_prime_3}, Val{:mu_prime_4}, Val{:mu_prime_5}, Val{:mu_prime_6}, Val{:mu_prime_9}, Val{:mu_prime_10}, Val{:mu_prime_11}, Val{:mu_prime_12}, Val{:mu_prime_15}, Val{:mu_prime_16}, Val{:mu_prime_17}, Val{:mu_prime_18}, Val{:mu_prime_19}, Val{:mu_prime_22}, Val{:mu_prime_23}, Val{:mu_prime_24}, Val{:mu_prime_25}})
#specs.distributions
ana_spec = specs.analyses[1]
a = HS3.make_analyses(ana_spec, specs)
f = unique!(reduce(vcat, a.free_parameters))
using ValueShapes
using Random
mmn = (; (i => Uniform(0, 1) for i in f)...)
prior = NamedTupleDist(mmn)
valus = rand(prior)
logdensityof(a.likelihood, valus)
k = NTuple{}()
typeof(k)
map(x -> (x),  a.free_parameters[i] for i in 1:11)
e = [merge(k, a.free_parameters[i]) for i in 1:11] 
println(keys(a.prior))
println((a.parameters_of_interest))
#@benchmark HS3.make_analyses(ana_spec, specs)
println(typeof(a.likelihood))
b = a.likelihood
using BAT
using Random
using Distributions, ValueShapes
b = rand(a.prior)
b = merge(b, (staterror= 0.,))

c = NamedTupleDist( a= Uniform(0,1), b= Uniform(2, 3))
keys(b)
posterior = PosteriorMeasure(a.likelihood, a.prior)
samples = bat_sample(posterior, MCMCSampling(mcalg = MetropolisHastings(), nsteps = 10^5, nchains = 4)).result
bat_findmode(a.likelihood)
#z
#b = merge(a.distributions, a.functions)
#HS3.topological_sort(b)
#
#####Distributions einlesen 
#data_spec = AbstractDistSpec(example_dict.distributions.data_dist)
#v = (param_lambda = 2,)
#model = make_markov_kernel(data_spec)
#model(v)
#
#prior_spec = AbstractDistSpec(example_dict.distributions.prior_dist)
#prior_model = make_markov_kernel(prior_spec)
#p = prior_model(v)
#
#example_dict[:priors]
#subst_prior(example_dict[:priors], (prior_dist = p,))