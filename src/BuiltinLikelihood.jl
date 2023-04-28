function generate_likelihood(dist, data::StatsBase.Histogram, params)
    #for  i in 1:length(bin_centers)
    #    bin_centers[i] = bin_centers[i] - 0.5
    #end    
    bin_edges = data.edges[1]
    bin_edges_left = bin_edges[1:end-1]
    bin_edges_right = bin_edges[2:end]
    bin_widths = bin_edges_right - bin_edges_left
    bin_centers = (bin_edges_right + bin_edges_left) / 2 
    lpdf = 0.0
    idx = eachindex(data.weights)
    for  i in 1:length(data.weights)
         #bin_centers[i] = round(bin_centers[i],  digits=0) 
         expected_counts = bin_widths[i] * Distributions.pdf(dist(params), bin_centers[i])
         #println(bin_centers)
         lpdf += Distributions.logpdf(Distributions.Poisson(expected_counts), data.weights[i])
    end
    #lpdf = Distributions.logpdf(dist(params), data.weights)
     return lpdf
end
export generate_likelihood

function generate_likelihood(dist, params)
## return likelihood for distributions with implied data, as used for aux_disitributions 
    if isa(dist.inner.a, HS3.DistributionSpec{:gaussian_dist})
        return logpdf(dist(params), 0)
    elseif isa(dist.inner.a, HS3.DistributionSpec{:poisson_dist})
        return logpdf(dist(params), 1)
    elseif isa(dist.inner.a, HS3.DistributionSpec{:lognormal_dist})
        return logpdf(dist(params), 1)
    else
        @warn "Specified auxiliary distribution is not associated with implied data"
    end
end


abstract type AbstractLikelihoodSpec end

Parameters.@with_kw struct LikelihoodSpec <: AbstractLikelihoodSpec
    data::AbstractArray{String}
    distributions::AbstractArray{String}
    aux_distributions::Union{AbstractArray{String}, Nothing} = nothing
end

function generate_likelihoodspec(dict::AbstractDict)
    content_dict = filter(entry -> entry[1] != :name, dict)
    @assert length(content_dict[:data]) == length(content_dict[:distributions]) "Likelihood $(dict[:name]) not defined correctly."
    LikelihoodSpec(; NamedTuple(content_dict)...)
end

function make_likelihood_specs(arr::AbstractArray)
    nt = NamedTuple()
    for element in arr
        nt = (; nt..., Symbol(element.name) => generate_likelihoodspec(element),) 
    end
    nt
end

export make_likelihood_specs

function make_likelihood(likelihood_spec::AbstractLikelihoodSpec, functional_specs::NamedTuple, data_specs::NamedTuple)
    generated_dist = []
    generated_data = []
    expected_counts = []
    sorted_specs = topological_sort(functional_specs)
    for i in 1:length(likelihood_spec.distributions)
        generated_dist = push!(generated_dist, make_functional(functional_specs[Symbol(likelihood_spec.distributions[i])], sorted_specs))
        generated_data = push!(generated_data, make_data(data_specs[Symbol(likelihood_spec.data[i])]))
    end
    println("was ", generated_dist, " was")
    DensityInterface.logfuncdensity(function (params)
        ll = 0.0
        for i in 1:length(likelihood_spec.distributions)
            ll += generate_likelihood(generated_dist[i], generated_data[i], params)
        end
        if likelihood_spec.aux_distributions !== nothing
            for aux in likelihood_spec.aux_distributions
                ll += generate_likelihood(functional_specs[Symbol(aux)], params)
            end
        end
        return ll
    end)
end

function make_likelihood(specs::NamedTuple, functional_specs::NamedTuple, data_specs::NamedTuple)
    nt = NamedTuple()
    for (k, v) in zip(keys(specs), specs)
        nt = merge(nt, (k => make_likelihood(v, functional_specs, data_specs),))
    end
    return nt
end

#function generate_likelihood(dist, data::StatsBase.Histogram)
#    DensityInterface.logfuncdensity(function (params)
#        expected = rand(Poisson(params))
#        lpdf = Distributions.logpdf(dist(params), data.weights[1])
#        idx = eachindex(data.weights)
#        for  i in 2:length(data.weights)
#             lpdf += Distributions.logpdf(dist(params), data.weights[i])
#        end
#        #lpdf = Distributions.logpdf(dist(params), data.weights)
#         return lpdf
#    end)
#end