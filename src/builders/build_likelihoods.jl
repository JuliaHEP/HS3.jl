# This file is a part of HS3.jl, licensed under the MIT License (MIT).
"""
    generate_likelihood(dist, data::StatsBase.Histogram)

Generate a likelihood function for a given distribution and histogram data in StatsBase.Hsitogram format.

# Arguments
- `dist`: The distribution.
- `data`: The histogram data.

# Returns
A likelihood function.
"""
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

function generate_likelihood(dist, data::UnbinnedData)  
    likelihood = function (params)
        lpdf = 0.0
        for i in 1:length(data.weights)
            lpdf += Distributions.logpdf(Distributions.Poisson(Distributions.pdf(dist(params), data.entries[i][1])), data.weights[i])
        end
        return lpdf
    end
    return params -> likelihood(params)
end

"""
    generate_likelihood(dist)

Generate a likelihood function based on the provided distribution.

# Arguments
- `dist`: A functional object. Usually these are defined as auxiliary disitrbutions.

# Returns
- A function that calculates the likelihood of the given distribution.

# Details
This function takes a function object `dist` and returns a new function that calculates the likelihood of the distribution. The generated likelihood function accepts a parameter vector `params` as input and returns the logarithm of the probability density function (PDF).
If the variables `x` are a vector, the function returns the log PDF values for each element of `x`. If `x` is a scalar value, the function returns the log PDF value for that single value.
"""
function generate_likelihood(dist)
    return params -> begin
        x = dist(params)[2]
        logpdf_x = length(x) > 1 ? [a for a in x] : x[1]
        Distributions.logpdf(dist(params)[1], logpdf_x) 
    end
end
    
 
"""
generate_likelihood(dist::HistfactPDF, data::StatsBase.Histogram)

Generate a likelihood function for a given HistfactPDF distribution and histogram data.

# Arguments
- `dist`: The HistfactPDF distribution.
- `data`: The histogram data.

# Returns
A likelihood function.
"""   
function generate_likelihood(dist::HistfactPDF, data::StatsBase.Histogram)  
    likelihood = LiteHF.pyhf_loglikelihoodof(dist.channel[1], convert.(Float64, data.weights))
    priors = LiteHF.pyhf_logpriorof(dist.prior)
    custom = Base.Fix1(_include_customs, dist)
    @info ((dist.prior))
    return (params  -> begin 
        params = custom(params)
        return likelihood([params[x] for x in dist.order]) + priors([params[x] for x in keys(dist.prior)])
        #return likelihood([params[x] for x in dist.order]) #+ priors([params[x] for x in keys(dist.prior)]) 
    end)
end
    

#@generated function _include_customs(order::Vector,customs::NamedTuple, params)
#    params_expr = :(params)
#    order_expr = :(order)
#    customs_expr = :(customs)
#    return :(($order_expr, $customs_expr, $params_expr) -> begin
#        for x in $order_expr
#            if haskey($customs_expr, x)
#                nt = (;)
#                for k in $customs_expr[x]
#                    nt = merge_namedtuples(nt, (k => $params_expr[k],))
#                end
#                $params_expr = merge_namedtuples($params_expr, (x => nt,))
#            end
#        end
#        return $params_expr
#    end)
#end

function _include_customs(dist, params)
        for x in dist.order
            if haskey(dist.customs, x)
                nt = (;)
                for k in dist.customs[x]
                    nt = merge(nt, (k => (params[k]),))
                end
                #if !haskey(params, x)
                    params = merge(params, (x => nt,))
                #end
            end
        end
    return params
end

"""
    make_likelihood(likelihood_spec::LikelihoodSpec, functional_specs::NamedTuple, data_specs::NamedTuple)

Create a likelihood function based on the given likelihood specification, functional specifications, and data specifications.

# Arguments
- `likelihood_spec`: The specification of the likelihood, including distributions, data and auxiliary distributions.
- `functional_specs`: The named tuple of functional specifications.
- `data_specs`: The named tuple of data specifications.

# Returns
A combined likelihood function.

"""
function make_likelihood(likelihood_spec::LikelihoodSpec, functional_specs::NamedTuple, data_specs::NamedTuple)
    generated_dist = []
    generated_data = []
    free_parameters = []
    sorted_specs = topological_sort(functional_specs)

    ### making distributions
    if likelihood_spec.distributions !== nothing
        for i in 1:length(likelihood_spec.distributions)
            if typeof(functional_specs[Symbol(likelihood_spec.distributions[i])]) == HistFactorySpec
                funct = make_histfact(functional_specs[Symbol(likelihood_spec.distributions[i])], Symbol(likelihood_spec.distributions[i]), functional_specs)
                #free_parameters = push!(free_parameters, collect(keys(funct.prior)))
                generated_dist = push!(generated_dist, funct)
            else
                generated_dist = push!(generated_dist, make_functional(functional_specs[Symbol(likelihood_spec.distributions[i])], sorted_specs))
            end
            generated_data = push!(generated_data, make_data(data_specs[Symbol(likelihood_spec.data[i])]))
        end
    end

    likelihood_functions = []
    for i in 1:length(generated_dist)
        likelihood_functions = push!(likelihood_functions, generate_likelihood(generated_dist[i], generated_data[i]))
    end
    #auxiliary disitributions 
    if likelihood_spec.aux_distributions !== nothing
        for aux in likelihood_spec.aux_distributions
            likelihood_functions = push!(likelihood_functions, generate_likelihood(make_functional(functional_specs[Symbol(aux)], sorted_specs)))
        end
    end
    return DensityInterface.logfuncdensity(
        function(params)
            ll = 0.0
            for ll_function in likelihood_functions
                ll += ll_function(params) 
            end
        return ll
    end)
end


"""
    make_likelihoods(specs::NamedTuple, functional_specs::NamedTuple, data_specs::NamedTuple)

Create multiple likelihood functions based on the given specifications, functional specifications, and data specifications.

# Arguments
- `specs`: The named tuple of likelihood specifications.
- `functional_specs`: The named tuple of functional specifications.
- `data_specs`: The named tuple of data specifications.

# Returns
A named tuple of likelihood functions, where the keys correspond to the keys in `specs` and the values are the corresponding likelihood functions.

"""
function make_likelihoods(specs::NamedTuple, functional_specs::NamedTuple, data_specs::NamedTuple)
    nt = NamedTuple()
    for (k, v) in zip(keys(specs), specs)
        nt = merge(nt, (k => make_likelihood(v, functional_specs, data_specs),))
    end
    return nt
end