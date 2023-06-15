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

"""
    generate_likelihood(dist)

Generate a likelihood function for a given distribution.

# Arguments
- `dist`: A distribution function that takes `params` as input and returns a tuple `(distribution, position vector)` representing the value and additional arguments for the distribution.

# Returns
A function that computes the logarithm of the probability density function (PDF) of the distribution.
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
    likelihood = LiteHF.pyhf_logjointof(dist.channel[1], data.weights, dist.prior)
    println(keys(dist.prior))
    return params -> likelihood([convert(Float64, params[x]) for x in keys(dist.prior)])
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
    #println("sorted: ", sorted_specs)
    ### making distributions
    if !isempty(sorted_specs)
        for i in 1:length(likelihood_spec.distributions)
            if typeof(functional_specs[Symbol(likelihood_spec.distributions[i])]) == HistFactorySpec
                funct = make_histfact(functional_specs[Symbol(likelihood_spec.distributions[i])], Symbol(likelihood_spec.distributions[i]))
                free_parameters = push!(free_parameters, collect(keys(funct.prior)))
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
    #auxiliary likelihoods 
    if likelihood_spec.aux_distributions !== nothing
        for aux in likelihood_spec.aux_distributions
            likelihood_functions = push!(likelihood_functions, generate_likelihood(make_functional(functional_specs[Symbol(aux)], sorted_specs)))
        end
    end
    return (likelihood = DensityInterface.logfuncdensity(
        function(params)
            ll = 0.0
            for ll_function in likelihood_functions
                ll += ll_function(params) 
            end
        return ll
    end), 
    free_parameters = unique!(free_parameters))
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