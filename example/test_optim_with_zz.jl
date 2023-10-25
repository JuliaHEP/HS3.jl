using HS3
using DensityInterface
using BenchmarkTools
using NamedTupleTools
using ValueShapes

#init 
dict = file_to_dict("/ceph/groups/e4/users/cburgard/public/forRobin/higgs/ouput/noES_ind/zz/125.json")
specs = HS3.generate_specs(dict)
analysis = HS3.make_analyses(specs.analyses[5], specs)

prior = HS3.generate_uniform_prior_from_domain(analysis[:parameter_domain])
constant_points = HS3.make_const_parameterpoints(specs.parameter_points.default_values)
all_points = HS3.make_parameterpoints(specs.parameter_points.default_values)

#remove unnecassary variables
prior = delete(NamedTuple(prior), (:obs_x_ATLAS_H_2e2mu_channel_2011_1250_1 , :obs_x_ATLAS_H_2e2mu_channel_2012_1250_1, :obs_x_ATLAS_H_2mu2e_channel_2011_1250_1, :obs_x_ATLAS_H_2mu2e_channel_2012_1250_1, :obs_x_ATLAS_H_4e_channel_2011_1250_1, :obs_x_ATLAS_H_4e_channel_2012_1250_1, :obs_x_ATLAS_H_4mu_channel_2011_1250_1, :obs_x_ATLAS_H_4mu_channel_2012_1250_1))
prior = merge(NamedTupleDist(prior), constant_points)
#prior = NamedTupleDist(prior)
totalndof(valshape(prior)) #should be 63
typeof(values(prior))
all_points
using BAT
posterior = PosteriorMeasure(analysis[:likelihood], (prior))
logdensityof(analysis[:likelihood], all_points)
logdensityof(posterior, all_points)
@btime logdensityof(analysis[:likelihood], all_points)
#init values
#x0 = collect(select(all_points, Tuple(setdiff(keys(all_points), keys(constant_points)))))
#x0 = Float64.(x0)
using Distributions
using ValueShapes
b = Normal(0, 0.036)
logpdf(b, 2) 
function reorder_nt(reference_nt, target_nt::NamedTuple)::NamedTuple
    return NamedTuple{keys(reference_nt)}(tuple([getfield((target_nt), k) for k in keys(reference_nt)]...))
end
all_points = merge(all_points, (mu =0.,))
all_points = reorder_nt(prior, all_points)

function unique_fields_values(nt1::NamedTuple, nt2::NamedTuple)
    unique_keys = setdiff(fieldnames(typeof(nt1)), fieldnames(typeof(nt2)))
    return NamedTuple{Tuple(unique_keys)}(tuple(map(k -> getproperty(nt1, k), unique_keys)...))
end
un = unique_fields_values(analysis[:parameter_domain], (constant_points))
un = delete(un,  (:obs_x_ATLAS_H_2e2mu_channel_2011_1250_1 , :obs_x_ATLAS_H_2e2mu_channel_2012_1250_1, :obs_x_ATLAS_H_2mu2e_channel_2011_1250_1, :obs_x_ATLAS_H_2mu2e_channel_2012_1250_1, :obs_x_ATLAS_H_4e_channel_2011_1250_1, :obs_x_ATLAS_H_4e_channel_2012_1250_1, :obs_x_ATLAS_H_4mu_channel_2011_1250_1, :obs_x_ATLAS_H_4mu_channel_2012_1250_1)) 

function make_namedtuple_from_file(filepath::String)
    pairs = Pair{Symbol, Float64}[]
    open(filepath, "r") do file
        for line in eachline(file)
            if occursin("=", line)
                parts = split(line, "=")
                var_name = strip(parts[1])
                var_value = parse(Float64, strip(parts[2]))
                push!(pairs, Symbol(var_name) => var_value)
            end
        end
    end
    return NamedTuple(pairs)
end

# Example usage
filepath = "/net/e4-nfs-home.e4.physik.tu-dortmund.de/home/rpelkner/Downloads/distribution_results.txt"
default = make_namedtuple_from_file(filepath)

ex = explore_values(un, analysis[:parameter_domain], analysis[:likelihood], default)

write_to_csv(ex, "zz_final.csv")


include 





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
    #println(f isa Function, typeof(x_init))
    #println(algorithm)
    #asdsaxc
    optim_result = my_optim_minimize(f, x_init, Optim.LBFGS())
    r_optim = Optim.MaximizationWrapper(optim_result)
    transformed_mode = Optim.minimizer(r_optim.res)
    result_mode = inv_trafo(transformed_mode)
#
#
#
    (result = result_mode, result_trafo = transformed_mode, trafo = trafo, #=trace_trafo = trace_trafo,=# info = r_optim)
    #sleep(10)
end

function my_optim_minimize(f, x_init, algorithm)
    @info "using local version of _optim_minimize for zeroth order optimizers!"
    algorithm = Optim.LBFGS(; m=100)
    opts = Optim.Options(store_trace = true, extended_trace=true, iterations = 4000, allow_f_increases=true, g_tol=1e-08 )
    Optim.optimize(f, x_init, algorithm, opts)
end

optimizer = OptimAlg(init = ExplicitInit([rand(prior)]), trafo=DoNotTransform())
mode1 = bat_findmode(posterior, optimizer)
keys(mode1)
mode1.info
mode1.result.mu
all_points
##check for minimum
logdensityof(posterior)(mode1.result)
logdensityof(posterior)(merge(mode1.result, (mu=0,)))
f = x -> logdensityof(posterior)(merge(mode1.result, (mu=x,)))
x = LinRange(0, .000000000001, 1000)
y = []
for val in x
    #print(" ", val)    
    push!(y, f(val))
end
y
using Plots 
plot(x, y)
#vline!([mode1.result.mu], label="Best Fit")



logdensityof(posterior)(rand(prior))


#init=ExplicitInit([default]),
using Profile

#speed tests 
analysis = HS3.make_analyses(specs.analyses[5], specs)
posterior = PosteriorMeasure(analysis[:likelihood], prior)

@btime logdensityof(posterior)(all_points)

using Profiler
@profview logdensityof(posterior)(all_points)
Profile.print(sortedby = :count)


@profview logdensityof(analysis[:likelihood])(all_points)
@btime logdensityof(analysis[:likelihood])(all_points)



@code_warntype logdensityof(posterior)(all_points)
using Distributed
@distributed for i in 1:10000
    println(i)
end

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
problem = LikelihoodProblem(loglik, x0, syms=names, f_kwargs=(adtype=Optimization.AutoFiniteDiff(),))

sh = varshape(posterior.prior)
f = (x, data) -> logdensityof(posterior)(sh(x))
f(x0)
typeof(collect(names))
problem = LikelihoodProblem(f, x0, syms=collect(names),  f_kwargs=(adtype=Optimization.AutoFiniteDiff(),))

# Maximum likelihood estimate (should be same result as bat_findmode)
sol = mle(problem, Optim.LBFGS(), ) # use Optim.NelderMead() to not use gradient information
for i in eachindex(sol.mle)
    println(names[i], " ", sol.mle[i])
end
# Profiling
lower_bounds = posterior.prior.bounds.vol.lo
upper_bounds = posterior.prior.bounds.vol.hi
resolutions = ones(Int, length(x0)) * 200 # one entry for each parameter

param_ranges = construct_profile_ranges(sol, lower_bounds, upper_bounds, resolutions)

prof = profile(problem, sol; conf_level=0.683, param_ranges, parallel=true)



prof
dump(prof)
prof.confidence_intervals[5][1]
collect(names)

namedtuple_entries = [(Symbol(name),(prof.confidence_intervals[key][1],  prof.confidence_intervals[key][2])) for (name, key) in zip(collect(names), sort(collect(keys(prof.confidence_intervals))))]
result = NamedTuple(namedtuple_entries)


using CSV
using JSON3
#JSON3.write("zz_profile_ll.txt", prof.confidence_intervals)
open("higgs_noES_zzCI.json", "w") do f
    JSON3.write(f, (result))
end


write("higgs_noES_zzCI.txt", result)
a=1

using ΒAT
samples = bat_sample(posterior)
bat_report(samples.result)
prior.alpha_pdf_gg_llll_2011
evm = EvaluatedMeasure(posterior, samples=samples.result)
write("report.txt", bat_report(evm))
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
get_smallest_interval_edges(samples.result, :mu, 0.683, atol=0.1)

using Plots
Plots.plot(samples.result)