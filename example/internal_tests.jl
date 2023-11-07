#NOTE: This is an INTERNAL testing file for developement purposes! 
#Use histfactory_hl.jl instead!

using HS3
using DensityInterface
using BenchmarkTools
using NamedTupleTools
using ValueShapes
using BAT
include("/net/e4-nfs-home.e4.physik.tu-dortmund.de/home/rpelkner/Documents/higgs_analyses/ram_sampler_v3.jl")
#init 
dict = file_to_dict("/net/e4-nfs-home.e4.physik.tu-dortmund.de/home/rpelkner/Downloads/ouput/noES_ind/zz/zz.json")
specs = HS3.generate_specs(dict)
analysis = HS3.make_analyses(specs.analyses[5], specs)



using Distributions
prior = HS3.generate_uniform_prior_from_domain(analysis[:parameter_domain])
var_names = [:ATLAS_EFF_SIG_4l, :ATLAS_E_CUTS_EFF, :ATLAS_E_EFF, :ATLAS_E_TRIG, :ATLAS_M_EFF, :ATLAS_M_TRIG, :ATLAS_SHAPE_BKG_4l, :BR_VV, :QCDscale_VH, :QCDscale_VV, :QCDscale_ggH, :QCDscale_ggVV, :QCDscale_qqH, :alpha_ATLAS_Norm_Z_2mu2e_2011_llll_2011, :alpha_ATLAS_Norm_Z_2mu2e_2012_llll_2012, :alpha_ATLAS_Norm_Z_4e_2011_llll_2011, :alpha_ATLAS_Norm_Z_4e_2012_llll_2012, :alpha_ATLAS_Norm_Zbb_2e2mu_2011_llll_2011, :alpha_ATLAS_Norm_Zbb_2e2mu_2012_llll_2012, :alpha_ATLAS_Norm_Zbb_4mu_2011_llll_2011, :alpha_ATLAS_Norm_Zbb_4mu_2012_llll_2012, :alpha_ATLAS_Norm_ttbar_2e2mu_2011_llll_2011, :alpha_ATLAS_Norm_ttbar_2e2mu_2012_llll_2012, :alpha_ATLAS_Norm_ttbar_4mu_2011_llll_2011, :alpha_ATLAS_Norm_ttbar_4mu_2012_llll_2012, :lumi2011, :lumi2012, :mu, :pdf_gg, :pdf_qqbar]
nt_vals= []
for i in var_names
    push!(nt_vals, Distributions.Uniform(first(analysis[:parameter_domain][i]), last(analysis[:parameter_domain][i])))
end
prior = NamedTuple{Tuple(var_names)}(Tuple(nt_vals))
prior = merge(all_points, prior)
#prior_dist = delete(NamedTuple(prior_dist), (:obs_x_ATLAS_H_2e2mu_channel_2011_1250_1 , :obs_x_ATLAS_H_2e2mu_channel_2012_1250_1, :obs_x_ATLAS_H_2mu2e_channel_2011_1250_1, :obs_x_ATLAS_H_2mu2e_channel_2012_1250_1, :obs_x_ATLAS_H_4e_channel_2011_1250_1, :obs_x_ATLAS_H_4e_channel_2012_1250_1, :obs_x_ATLAS_H_4mu_channel_2011_1250_1, :obs_x_ATLAS_H_4mu_channel_2012_1250_1))
prior = NamedTupleDist(prior)
constant_points = HS3.make_const_parameterpoints(specs.parameter_points.default_values)
all_points = HS3.make_parameterpoints(specs.parameter_points.default_values)

#remove unnecassary variables
prior = delete(NamedTuple(prior), (:obs_x_ATLAS_H_2e2mu_channel_2011_1250_1 , :obs_x_ATLAS_H_2e2mu_channel_2012_1250_1, :obs_x_ATLAS_H_2mu2e_channel_2011_1250_1, :obs_x_ATLAS_H_2mu2e_channel_2012_1250_1, :obs_x_ATLAS_H_4e_channel_2011_1250_1, :obs_x_ATLAS_H_4e_channel_2012_1250_1, :obs_x_ATLAS_H_4mu_channel_2011_1250_1, :obs_x_ATLAS_H_4mu_channel_2012_1250_1))

#set constant values constant
prior = merge((prior), constant_points)
prior = NamedTupleDist(prior)
prior = NamedTupleDist(prior)
#test performance
params = rand(prior)
logdensityof(analysis[:likelihood], params)



#alternatively I use this to build custom simpler HistFactory models
dict = file_to_dict("example/simpleHistFact.json")
specs = HS3.generate_specs(dict)
analysis = HS3.make_analyses(specs.analyses[1], specs)
logdensityof(analysis[:likelihood], (mu=1., syst1=0., syst2=0.0, syst3 = 0., lumi = 1., ))




#sort the namedtuple
function reorder_nt(reference_nt, target_nt::NamedTuple)::NamedTuple
    return NamedTuple{keys(reference_nt)}(tuple([getfield((target_nt), k) for k in keys(reference_nt)]...))
end

all_points = merge(all_points, (mu =0.,))
all_points = reorder_nt(prior, all_points)

totalndof(valshape(prior)) #should be 63


posterior = PosteriorMeasure(analysis[:likelihood], (prior))

#@elapsed samples = bat_sample(posterior, MCMCSampling(mcalg = MetropolisHastings(), nsteps = 10^1, nchains = 4)).result
#
#@elapsed samples = bat_sample(posterior, MCMCSampling(mcalg = MetropolisHastings(), nsteps = 10^4, nchains = 4)).result
using Distributions
algorithm = RAMSampler(nchains=2, nsteps=10^6, nburnin=10^5)
@elapsed samples = bat_sample(posterior, algorithm)
println(bat_report(samples.result))

using StatsBase
function get_smallest_interval_edges(
    samples::DensitySampleVector,
    key::Union{Symbol,Real},
    p::Real;
    bins=200,
    atol=0.0)

    marg = BAT.MarginalDist(samples, key, bins=bins)

    marghist = convert(Histogram, marg.dist)
    hist_p = BAT.get_smallest_intervals(marghist, [p])[1][1]

    lower, upper = BAT.get_interval_edges(hist_p, atol=atol)
    return (lower=lower, upper=upper)
end
#varshape(prior_dist)

get_smallest_interval_edges(samples.result, :mu, 0.9999998 , bins=400, atol=0.3)

for i in var_names
    println(i, " " , get_smallest_interval_edges(samples.result, i, 0.682689492 , bins=200, atol=0.0))
end

#overwrite bat_findmode
import BAT
import BAT: bat_findmode
using Optim
using BAT: AnyMeasureOrDensity, AbstractMeasureOrDensity
using BAT: get_context, get_adselector, _NoADSelected
using BAT: bat_initval, transform_and_unshape, apply_trafo_to_init
using BAT: negative, inverse, fchain
using DensityInterface, ChangesOfVariables, InverseFunctions, FunctionChains

function BAT.bat_findmode_impl(target::AnyMeasureOrDensity, algorithm::OptimAlg, context::BATContext)
   @info "Using my own bat_findmode_impl"
    transformed_density, trafo = transform_and_unshape(algorithm.trafo, target, context)
    inv_trafo = inverse(trafo)

    initalg = apply_trafo_to_init(trafo, algorithm.init)
    x_init = collect(bat_initval(transformed_density, initalg, context).result)
    # Maximize density of original target, but run in transformed space, don't apply LADJ:
    f = fchain(inv_trafo, logdensityof(target), -)
    optim_result = my_optim_minimize(f, x_init, Optim.LBFGS())
    r_optim = Optim.MaximizationWrapper(optim_result)
    transformed_mode = Optim.minimizer(r_optim.res)
    result_mode = inv_trafo(transformed_mode)
    (result = result_mode, result_trafo = transformed_mode, trafo = trafo, #=trace_trafo = trace_trafo,=# info = r_optim)
end

function my_optim_minimize(f, x_init, algorithm)
    @info "using local version of _optim_minimize for zeroth order optimizers!"
    algorithm = Optim.LBFGS(; m=10)
    opts = Optim.Options(store_trace = true, extended_trace=true, iterations = 1000, allow_f_increases=true, g_tol=1e-08 )
    Optim.optimize(f, x_init, algorithm, opts)
end

optimizer = OptimAlg(init = ExplicitInit([rand(prior)]), trafo=PriorToGaussian())
@elapsed mode1 = bat_findmode(posterior, optimizer)
mode1.info
mode1.result
 

using ProfileLikelihood 
#docu: https://danielvandh.github.io/ProfileLikelihood.jl/dev/


var_names = filter(x -> x in var_names, (fieldnames(typeof(NamedTuple(prior)))))

# make the BAT posterior compatible with ProfileLikelihood.jl
llikelihood = logdensityof(unshaped(posterior.likelihood))
loglik(x, data) = llikelihood(x)
names = collect(keys(varshape(posterior)))
#confidence_interval_method=:extrema
using Optimization
using OptimizationOptimJL
import ForwardDiff
using AutoDiffOperators
x0 = collect(select(rand(prior), Tuple(setdiff(keys(all_points), keys(constant_points)))))
x0 = Float64.(x0)
names = keys(select(rand(prior), Tuple(setdiff(keys(all_points), keys(constant_points)))))
#x0 = [0., 0., 0.] # starting value for the minimization
# create ProfileLikelihood problem
sh = varshape(posterior.prior)
f = (x, data) -> logdensityof(posterior)(sh(x))
problem = LikelihoodProblem(f, x0; syms=collect(names),  f_kwargs=(adtype=Optimization.AutoFiniteDiff(),),)

problem.syms
# Maximum likelihood estimate (should be same result as bat_findmode)
#sh = varshape(posterior.prior)
println(totalndof(varshape(posterior.prior)))

sol = mle(problem, Optim.LBFGS(; m=100),) # use Optim.NelderMead() to not use gradient information
for i in eachindex(sol.mle)
    println(names[i], " ", sol.mle[i])
end
# Profiling
lower_bounds = posterior.prior.bounds.vol.lo
upper_bounds = posterior.prior.bounds.vol.hi
resolutions = ones(Int, length(x0)) * 50 # one entry for each parameter
 
param_ranges = construct_profile_ranges(sol, lower_bounds, upper_bounds, resolutions)
σ = 0.682689492 
prof = profile(problem, sol; conf_level=σ, param_ranges, parallel=true, normalise = false)
#prof = profile(problem, sol; conf_level=σ, param_ranges, parallel=true)
typeof(sol)
typeof(problem)