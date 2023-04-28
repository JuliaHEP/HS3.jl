abstract type AbstractDomainsSpec <: AbstractHS3Spec end

function generate_domainspecs(::Val{:product_domain}, arr::AbstractArray)
    nt = NamedTuple()
    for element in arr
        temp = NamedTuple(filter(entry -> entry[1] != :name, element))
        nt = merge(nt, (Symbol(element.name) => DomainSpec(; temp...),))
    end
    return nt
end

Parameters.@with_kw struct DomainSpec <: AbstractDomainsSpec 
    max::Number
    min::Number 
end

function make_domainspecs(domain_array::AbstractArray)
    specs = NamedTuple()
    for domain in domain_array
        specs = merge(specs, (Symbol(domain.name) => generate_domainspecs(Val(Symbol(domain.type)), domain.axes),))
    end
    specs
end

function make_domain(domain_specs::NamedTuple)
    domain = NamedTuple()
    for (k, v) in zip(keys(domain_specs), domain_specs)
        print(k, v)
        domain = merge(domain, (Symbol(k) => (v.min:(v.max-v.min):v.max),))
    end
    domain
end

function generate_prior_from_domain(domain_specs::NamedTuple)
    prior = NamedTuple()
    for (k, v) in zip(keys(domain_specs), domain_specs)
        prior = merge(prior, (Symbol(k) => (Distributions.Uniform(v.min, v.max)),))
    end
    ValueShapes.NamedTupleDist(prior)
end
