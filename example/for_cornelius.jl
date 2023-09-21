using HS3
using DensityInterface
using BenchmarkTools
using NamedTupleTools
using ValueShapes
using BAT

#init 
dict = file_to_dict("/ceph/groups/e4/users/cburgard/public/forRobin/higgs/ouput/noES_ind/zz/125.json")
specs = HS3.generate_specs(dict)
analysis = HS3.make_analyses(specs.analyses[5], specs)


prior = HS3.generate_uniform_prior_from_domain(analysis[:parameter_domain])
constant_points = HS3.make_const_parameterpoints(specs.parameter_points.default_values)
all_points = HS3.make_parameterpoints(specs.parameter_points.default_values)

#remove unnecassary variables
prior = delete(NamedTuple(prior), (:obs_x_ATLAS_H_2e2mu_channel_2011_1250_1 , :obs_x_ATLAS_H_2e2mu_channel_2012_1250_1, :obs_x_ATLAS_H_2mu2e_channel_2011_1250_1, :obs_x_ATLAS_H_2mu2e_channel_2012_1250_1, :obs_x_ATLAS_H_4e_channel_2011_1250_1, :obs_x_ATLAS_H_4e_channel_2012_1250_1, :obs_x_ATLAS_H_4mu_channel_2011_1250_1, :obs_x_ATLAS_H_4mu_channel_2012_1250_1))

#set constant values constant
prior = merge(NamedTupleDist(prior), constant_points)
#prior = NamedTupleDist(prior)

#sort the namedtuple
function reorder_nt(reference_nt, target_nt::NamedTuple)::NamedTuple
    return NamedTuple{keys(reference_nt)}(tuple([getfield((target_nt), k) for k in keys(reference_nt)]...))
end

all_points = merge(all_points, (mu =0.,))
all_points = reorder_nt(prior, all_points)

totalndof(valshape(prior)) #should be 63


posterior = PosteriorMeasure(analysis[:likelihood], (prior))




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

optimizer = OptimAlg(init = ExplicitInit([rand(prior)]), trafo=DoNotTransform())
mode1 = bat_findmode(posterior, optimizer)
mode1.info
mode1.result.mu
 

using ProfileLikelihood 
#docu: https://danielvandh.github.io/ProfileLikelihood.jl/dev/

# make the BAT posterior compatible with ProfileLikelihood.jl
llikelihood = logdensityof(unshaped(posterior.likelihood))
loglik(x, data) = llikelihood(x)
names = collect(keys(varshape(posterior)))

using Optimization
using OptimizationOptimJL
x0 = collect(select(rand(prior), Tuple(setdiff(keys(all_points), keys(constant_points)))))
x0 = Float64.(x0)
names = keys(select(rand(prior), Tuple(setdiff(keys(all_points), keys(constant_points)))))
#x0 = [0., 0., 0.] # starting value for the minimization
# create ProfileLikelihood problem
sh = varshape(posterior.prior)
f = (x, data) -> logdensityof(posterior)(sh(x))
problem = LikelihoodProblem(f, x0, syms=collect(names),  f_kwargs=(adtype=Optimization.AutoFiniteDiff(),))

# Maximum likelihood estimate (should be same result as bat_findmode)
sol = mle(problem, Optim.LBFGS(), ) # use Optim.NelderMead() to not use gradient information
for i in eachindex(sol.mle)
    println(names[i], " ", sol.mle[i])
end
# Profiling
lower_bounds = posterior.prior.bounds.vol.lo
upper_bounds = posterior.prior.bounds.vol.hi
resolutions = ones(Int, length(x0)) * 5 # one entry for each parameter

param_ranges = construct_profile_ranges(sol, lower_bounds, upper_bounds, resolutions)

prof = profile(problem, sol; conf_level=0.683, param_ranges, parallel=true)