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
    println(bin_edges)
    bin_edges_left = bin_edges[1:end-1]
    bin_edges_right = bin_edges[2:end]
    bin_widths = bin_edges_right - bin_edges_left
    bin_centers = (bin_edges_right + bin_edges_left) / 2
    @info bin_widths bin_centers 
    likelihood = function (params)
        lpdf = 0.0
        for i in 1:length(data.weights)
            #@info Distributions.pdf(Distributions.truncated(dist(params), first(bin_edges), last(bin_edges)), bin_centers[i])
            #@info dist(params)  Distributions.logpdf(ContinousPoisson(bin_widths[i] * Distributions.pdf(Distributions.truncated(dist(params), first(bin_edges), last(bin_edges)), bin_centers[i])), data.weights[i]) 
            @info dist(params) bin_centers[i] first(bin_edges) last(bin_edges) 
            @info collect(bin_widths)
            lpdf += Distributions.logpdf(ContinousPoisson(bin_widths[i] * Distributions.pdf(Distributions.truncated(dist(params), first(bin_edges), last(bin_edges)), bin_centers[i])), data.weights[i]) 
        end
        return lpdf
    end
    return params -> likelihood(params) #+ 867.1179624924561
end

function generate_likelihood(dist, data::UnbinnedData)  
    #@info first(data.axes[1].range), last(data.axes[1].range)
    likelihood = function (params)
        range_low = first(data.axes[1].range) 
        range_high = last(data.axes[1].range)
        norm = log(_norm(dist(params), range_low, range_high))
        lpdf = 0.0
        for i in 1:length(data.entries)
            lpdf += (_normalize_logpdf(dist(params), range_low, range_high, data.entries[i][1]) - (norm)) * data.weights[i]
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
    likelihood = pyhf_loglikelihoodof(dist.channel, Float64.(data.weights))
    order = collect(dist.order)
    function ll(params::NamedTuple)::Float64  
        return (likelihood( [params[x] for x in order]))
    end
    return ll
end

function generate_constraints(constraint_dists::NamedTuple)
    constraints = pyhf_logpriorof(constraint_dists)

    constraint_names = collect(keys(constraint_dists))
    function ll(params::NamedTuple)::Float64  
        constraints(([params[x] for x in constraint_names]))
    end
    return ll
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
function make_likelihood(likelihood_spec::LikelihoodSpec, functional_specs::NamedTuple, data_specs::NamedTuple, domain::NamedTuple)
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
            generated_data = push!(generated_data, make_data(data_specs[Symbol(likelihood_spec.data[i])], domain))
        end
    end

    likelihood_functions = Function[]
    constraints = (;)
    for i in 1:length(generated_dist)
        likelihood_functions = push!(likelihood_functions, generate_likelihood(generated_dist[i], generated_data[i]))
        if typeof(generated_dist[i]) == HistfactPDF
            constraints = merge(constraints, generated_dist[i].prior)
        end
    end
    if constraints != (;)
        likelihood_functions = push!(likelihood_functions, generate_constraints(constraints))
    end
    #auxiliary disitributions 
    if likelihood_spec.aux_distributions !== nothing
        for aux in likelihood_spec.aux_distributions
            likelihood_functions = push!(likelihood_functions, generate_likelihood(make_functional(functional_specs[Symbol(aux)], sorted_specs)))
        end
    end
    #chunks = Iterators.partition(spec_dict, length(keys(spec_dict)) ÷ Threads.nthreads()+1)
    #tasks = map(chunks) do chunk
    #    Threads.@spawn _gen_specs(chunk)
    #end
    #chunk_sums = fetch.(tasks)
    #return merge(chunk_sums...)
    #if length(likelihood_functions) < 5
    #function log_density(params::NamedTuple)
    #    # Pre-allocate the results
    #    n = length(likelihood_functions)
    #    results = Vector{Float64}(undef, n)
    #    
    #    # Compute the log-likelihoods in-place
    #    for (i, ll_function) in enumerate(likelihood_functions)
    #        results[i] = ll_function(params)
    #    end
    #
    #    # Sum the results
    #    return sum(results)
    #end
    function log_density(params::NamedTuple; likelihood_functions::Vector{Function}=likelihood_functions)
        return sum(ll_function(params) for ll_function in likelihood_functions)
    end
    return DensityInterface.logfuncdensity(log_density)
    #else 
    #calc = Base.Fix1(compute_density, likelihood_functions)
    #return DensityInterface.logfuncdensity(
    #    params -> calc(params)
    #)
    #end
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