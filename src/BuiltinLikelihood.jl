function generate_likelihood(dist, data::StatsBase.Histogram)  
    bin_edges = data.edges[1]
    bin_edges_left = bin_edges[1:end-1]
    bin_edges_right = bin_edges[2:end]
    bin_widths = bin_edges_right - bin_edges_left
    bin_centers = (bin_edges_right + bin_edges_left) / 2 
    #idx = eachindex(data.weights)
 
    ll = function (params)
        lpdf = 0.0
        for i in 1:length(data.weights)
            lpdf += Distributions.logpdf(Distributions.Poisson(bin_widths[i] * Distributions.pdf(dist(params), bin_centers[i])), data.weights[i])
        end
        return lpdf
    end
    return params -> ll(params)
end
export generate_likelihood

function generate_likelihood(dist)
## return likelihood for distributions with implied data, as used for aux_disitributions 
    if isa(dist.inner.a, HS3.DistributionSpec{:gaussian_dist})
        return begin params -> logpdf(dist(params), 0) end
    elseif isa(dist.inner.a, HS3.DistributionSpec{:poisson_dist})
        return  begin params -> logpdf(dist(params), 1) end
    elseif isa(dist.inner.a, HS3.DistributionSpec{:lognormal_dist})
        return begin params -> logpdf(dist(params), 1) end
    else
        @warn "Specified auxiliary distribution $dist is not associated with implied data"
    end
end

function generate_likelihood(dist::HistfactPDF, data::StatsBase.Histogram)  
    lpdf = LiteHF.pyhf_logjointof(dist.channel[1], data.weights, dist.prior)
    return  params -> lpdf([convert(Float64, params[x]) for x in keys(dist.prior)])
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
    #functional_specs_without_histfact = NamedTupleTools.namedtuple(filter((x) -> (typeof(x[2]) != HistFactorySpec), NamedTupleTools.convert(Dict, functional_specs)))
    
    #histfacts = NamedTupleTools.namedtuple(filter((x) -> (typeof(x[2]) == HistFactorySpec), NamedTupleTools.convert(Dict, functional_specs)))

    #println(functional_specs_without_histfact)

    ### sorting
    sorted_specs = topological_sort(functional_specs)
    println("sorted: ", sorted_specs)
    ### making disitrbtions
    if !isempty(sorted_specs)
        for i in 1:length(likelihood_spec.distributions)
            if typeof(functional_specs[Symbol(likelihood_spec.distributions[i])]) == HistFactorySpec
                generated_dist = push!(generated_dist, make_histfact(functional_specs[Symbol(likelihood_spec.distributions[i])], Symbol(likelihood_spec.distributions[i])))
            else
                generated_dist = push!(generated_dist, make_functional(functional_specs[Symbol(likelihood_spec.distributions[i])], sorted_specs))
            end
            generated_data = push!(generated_data, make_data(data_specs[Symbol(likelihood_spec.data[i])]))
        end
    end
    #generated_hist = make_histfact(functional_specs[1], keys(functional_specs)[1])
    #println("Keys: ", keys(generated_hist.prior))
    #generate the individual likelihood functions aka distribution-data pairs
    likelihood_functions = []
    for i in 1:length(generated_dist)
        likelihood_functions = push!(likelihood_functions, generate_likelihood(generated_dist[i], generated_data[i]))
    end
    #histfact_ll = generate_likelihood(generated_hist, generated_data[1])
    #return params -> begin BAT.LogDVal(histfact_ll(params)) end
    return DensityInterface.logfuncdensity(
        function(params)
            ll = 0.0
            #ll = (likelihood_functions[i](params) for i in 1:length(likelihood_spec.distributions))
            for i in 1:length(likelihood_spec.distributions)
                ll += likelihood_functions[i](params) 
            end
            #generate_likelihood(generated_hist, generated_data[1]) #for i in 1:length(likelihood_spec.distributions) 
            #if likelihood_spec.aux_distributions !== nothing
            #    for aux in likelihood_spec.aux_distributions
            #        ll += generate_likelihood(functional_specs[Symbol(aux)])
            #    end
            #end
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