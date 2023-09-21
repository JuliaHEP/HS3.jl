using HS3
using DensityInterface
using BenchmarkTools
dict = file_to_dict("./example/test_interp0d.json")
specs = HS3.generate_specs(dict)

func = HS3.make_functional(specs.distributions.model_channel1, HS3.topological_sort(merge(specs.distributions, specs.functions)), 0)
using Distributions
func((mu=1, a= 2))
analysis = HS3.make_analyses(specs.analyses[1], specs)
a = logdensityof(analysis.likelihood, (mu =0.2, a=2,))
HS3.topological_sort(merge(specs.distributions, specs.functions))
using DensityInterface
a = logdensityof(analysis.likelihood)
specs.functions.mu_interp

a = (b =1, c=2)
keys(a)
a[:b]

using Distributions
a = MixtureModel([Exponential(2), Normal(0, 1)], [0.5, 0.5])
pdf(a, 2)


using HS3
using IMinuit

using DensityInterface
using BenchmarkTools
dict = file_to_dict("/net/e4-nfs-home.e4.physik.tu-dortmund.de/home/rpelkner/Documents/125.json")
specs = HS3.generate_specs(dict)
specs.likelihoods[5].distributions
specs.distributions.model_ATLAS_H_4mu_channel_2011_1250_1_llll_2011_sep_edit.samples.ATLAS_Bkg_gg_ZZ_ATLAS_H_4mu_channel_2011_1250_1_overallSyst_x_HistSyst_llll_2011_model.modifiers#.pdf_gg
using Random, ValueShapes


using Parsers

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
named_tuple = make_namedtuple_from_file(filepath)

using Minuit
using Distributions
n = truncated(Normal(1, 0.018), 0, 10)
pdf(n, 1)
func = HS3.make_functional(specs.distributions.model_channel1, HS3.topological_sort(merge(specs.distributions, specs.functions)), 0)
using Distributions
func((mu=1, a= 2))
analysis = HS3.make_analyses(specs.analyses[5], specs)
analysis[:likelihood]
keys(analysis)
named_tuple.pdf_gg
prior = HS3.generate_uniform_prior_from_domain(analysis[:parameter_domain])
c = HS3.make_const_parameterpoints(specs.parameter_points.default_values)
p = HS3.make_parameterpoints(specs.parameter_points.default_values)
setdiff(names, keys(p))
length(p)
using ProfileLikelihood 
using BAT
using DensityInterface, ValueShapes
#docu: https://danielvandh.github.io/ProfileLikelihood.jl/dev/
posterior = PosteriorMeasure(analysis[:likelihood], (prior))
# make the BAT posterior compatible with ProfileLikelihood.jl
llikelihood = logdensityof(unshaped(posterior))
l(x) =  logdensityof(unshaped(posterior), x)
loglik(x, data) = llikelihood(x)
names = collect(keys(varshape(posterior)))
p = merge(p, (obs_x_ATLAS_H_2e2mu_channel_2011_1250_1 = 999.75,
obs_x_ATLAS_H_2e2mu_channel_2012_1250_1 = 999.75,
obs_x_ATLAS_H_2mu2e_channel_2011_1250_1 = 999.75,
obs_x_ATLAS_H_2mu2e_channel_2012_1250_1 = 999.75,
obs_x_ATLAS_H_4e_channel_2011_1250_1 = 999.75,
obs_x_ATLAS_H_4e_channel_2012_1250_1 = 999.75,
obs_x_ATLAS_H_4mu_channel_2011_1250_1 = 999.75,
obs_x_ATLAS_H_4mu_channel_2012_1250_1 = 999.75,))
using Optimization
using OptimizationOptimJL
prior = merge(prior, c)
prior = delete(NamedTuple(prior), (:obs_x_ATLAS_H_2e2mu_channel_2011_1250_1 ,
:obs_x_ATLAS_H_2e2mu_channel_2012_1250_1,
:obs_x_ATLAS_H_2mu2e_channel_2011_1250_1,
:obs_x_ATLAS_H_2mu2e_channel_2012_1250_1,
:obs_x_ATLAS_H_4e_channel_2011_1250_1,
:obs_x_ATLAS_H_4e_channel_2012_1250_1,
:obs_x_ATLAS_H_4mu_channel_2011_1250_1,
:obs_x_ATLAS_H_4mu_channel_2012_1250_1))
prior =NamedTupleDist(prior)
x0 = collect()# starting value for the minimization
using NamedTupleTools
x0 = collect(select(p, Tuple(setdiff(keys(p), keys(c)))))
x0 = Float64.(x0)
totalndof(valshape(prior))
# create ProfileLikelihood problem
sh = varshape(posterior.prior)
f = x -> logdensityof(posterior)(sh(x))
problem = LikelihoodProblem(loglik, x0, syms=names, )
f(x0)
prior
logdensityof(posterior, p2) - logdensityof(posterior.likelihood, p2) 
p2 = rand(prior)
# Maximum likelihood estimate (should be same result as bat_findmode)
using Plots

sol = mle(problem, Optim.LBFGS(), )
sol = optimize(f, x0, NelderMead())
sol = optimize(f, x0, LBFGS())

sh(x0) == p
p = reorder_nt(prior, p)

function reorder_nt(reference_nt, target_nt::NamedTuple)::NamedTuple
    return NamedTuple{keys(reference_nt)}(tuple([getfield(target_nt, k) for k in keys(reference_nt)]...))
end

import BAT: _optim_minimize
function _optim_minimize(f::Function, x_init::AbstractArray{<:Real}, algorithm::Optim.ZerothOrderOptimizer, ::BATContext)
    @info "using local version of _optim_minimize for zeroth order optimizers!"
    opts = Optim.Options(store_trace = true, extended_trace=true)
    Optim.optimize(f, x_init, algorithm, opts)
end
_optim_minimize(f, x0, Optim.NelderMead(), set_batcontext())
bat_findmode(posterior)
a = logdensityof(analysis[:likelihood], (named_tuple) )
named_tuple
function unique_fields_values(nt1::NamedTuple, nt2::NamedTuple)
    unique_keys = setdiff(fieldnames(typeof(nt1)), fieldnames(typeof(nt2)))
    return NamedTuple{Tuple(unique_keys)}(tuple(map(k -> getproperty(nt1, k), unique_keys)...))
end


named_tuple1 = (a=1, b=2, c=3)
named_tuple2 = (c=4, d=5, e=6)
un = unique_fields_values(analysis[:parameter_domain], (c))
 
specs.parameter_points.default_values
nt = (pdf_gg = 2, )
value_nt = explore_values(un, analysis[:parameter_domain], analysis[:likelihood], named_tuple)
value_nt.ATLAS_E_EFF
using CSV
function write_to_csv(nt::NamedTuple, filename::String)
    # Create an array to collect rows for the CSV
    rows = []
    
    # Process each key-value pair in the NamedTuple
    for (key, values) in pairs(nt)
        for (param_value, func_output) in values
            push!(rows, (param_name=key, param_value=param_value, func_output=func_output))
        end
    end
    
    # Write the data to a CSV file
    CSV.write(filename, rows)
end

write_to_csv(value_nt, "output_zz.csv")

function explore_values(nt::NamedTuple, ranges::NamedTuple, ll, default::NamedTuple)
    # Initialize a Dict for collecting the results
    results_dict = Dict{Symbol, Vector{Tuple{Float64, Float64}}}()
    
    # For each key in the namedtuple
    for key in keys(nt)
        # If the key has a specified range
        if haskey(ranges, key)
            # A list to store tuples of (param_value, func_output) for this key
            key_results = []
            
            # Loop over the range
            for val in ranges[key]
                # Modify the namedtuple
                modified_nt = merge(default, (key => val,))
                
                # Compute the result
                result = logdensityof(ll, modified_nt) 
                
                # Store the parameter value and the result
                key_results = push!(key_results, (val, result))
            end
            
            # Save the results for this key in the results_dict
            results_dict[key] = key_results
        end
    end
    
    # Convert the Dict back to a NamedTuple for the final result
    return NamedTuple(results_dict)
end



length(c)
#length(rand(prior))
println(rand(prior))
prior = merge(prior, c)
setdiff(keys(prior), keys(c))
prior
using DensityInterface
using BAT
using ValueShapes
using Optim
import ForwardDiff
using AutoDiffOperators

typeof(rand(alls))
optimizer = OptimAlg(optalg=Optim.NelderMead(), trafo=PriorToGaussian(), init=ExplicitInit([default]), options=Optim.Options(
iterations = 3000))
global_mode = bat_findmode(posterior, optimizer)
fieldnames(global_mode)
(global_mode.info)






typeof(default) 
typeof(alls)
2+2
length(keys(alls)) == length(keys(default))
typeof(collect(values(default)))
alls = merge(prior, c)

alls
set_batcontext(ad = ADModule(:ForwardDiff), autodiff = :forward)
default = HS3.make_parameterpoints(specs.parameter_points.default_values)
collect(values(default))
setdiff(keys(alls), keys(default))
using NamedTupleTools
((alls))
default = merge(default, (obs_x_ATLAS_H_2e2mu_channel_2011_1250_1 = 999.75,
obs_x_ATLAS_H_2e2mu_channel_2012_1250_1 = 999.75,
obs_x_ATLAS_H_2mu2e_channel_2011_1250_1 = 999.75,
obs_x_ATLAS_H_2mu2e_channel_2012_1250_1 = 999.75,
obs_x_ATLAS_H_4e_channel_2011_1250_1 = 999.75,
obs_x_ATLAS_H_4e_channel_2012_1250_1 = 999.75,
obs_x_ATLAS_H_4mu_channel_2011_1250_1 = 999.75,
obs_x_ATLAS_H_4mu_channel_2012_1250_1 = 999.75,)) 
prior = merge(prior, (obs_x_ATLAS_H_2e2mu_channel_2011_1250_1 = 999.75,
obs_x_ATLAS_H_2e2mu_channel_2012_1250_1 = 999.75,
obs_x_ATLAS_H_2mu2e_channel_2011_1250_1 = 999.75,
obs_x_ATLAS_H_2mu2e_channel_2012_1250_1 = 999.75,
obs_x_ATLAS_H_4e_channel_2011_1250_1 = 999.75,
obs_x_ATLAS_H_4e_channel_2012_1250_1 = 999.75,
obs_x_ATLAS_H_4mu_channel_2011_1250_1 = 999.75,
obs_x_ATLAS_H_4mu_channel_2012_1250_1 = 999.75,)) 



default = reorder_nt(prior, default)
default.pdf_qqbar_ACCEPT
alls
#bat_mode_result = bat_findmode(posterior, OptimAlg(optalg = Optim.ParticleSwarm()))
#bat_mode_result2 = bat_findmode(posterior, OptimAlg(optalg = Optim.NelderMead()))
s = bat_sample(posterior)
named_tuple.obs_x_ATLAS_H_4mu_channel_2012_1250_1
HS3.topological_sort(merge(specs.distributions, specs.functions))
using DensityInterface
specs.functions.mu_interp

a = (b =1, c=2)
keys(a)
a[:b]

using Distributions
a = MixtureModel([Exponential(2), Normal(0, 1)], [0.5, 0.5])
pdf(a, 2)