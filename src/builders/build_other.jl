# This file is a part of HS3.jl, licensed under the MIT License (MIT).
"""
    make_domain_ranges(domain_specs::NamedTuple)

Create domain ranges from the given domain specifications.

# Arguments
- `domain_specs::NamedTuple`: The domain specifications as a named tuple.

# Returns
A named tuple representing the domain ranges.

# Example
```julia
domain_specs = (x = ProductDomainSpec(min = 0, max = 1), y = ProductDomainSpec(min = -1, max = 1))
domain_ranges = make_domain_ranges(domain_specs)
"""
function make_domain_ranges(domain_specs::NamedTuple)
    domain = NamedTuple()
    for (k, v) in zip(keys(domain_specs), domain_specs)
        domain = merge(domain, (Symbol(k) => (v.min:(v.max-v.min):v.max),))
    end
    domain
end


"""
    generate_uniform_prior_from_domain(domain_specs::NamedTuple)

Generate a uniform prior distribution from the given domain specifications.
*This function only exist for testing convenience. It might be removed in the future.*

# Arguments
- `domain_specs::NamedTuple`: The domain specifications as a named tuple.

# Returns
A `ValueShapes.NamedTupleDist` object representing the uniform prior distribution.

# Example
```julia
domain_specs = (x = (min = 0, max = 1), y = (min = -1, max = 1))
prior = generate_uniform_prior_from_domain(domain_specs)

"""
function generate_uniform_prior_from_domain(domain_specs::NamedTuple)
    prior = NamedTuple()
    for (k, v) in zip(keys(domain_specs), domain_specs)
        prior = merge(prior, (Symbol(k) => (Distributions.Uniform(v.min, v.max)),))
    end
    ValueShapes.NamedTupleDist(prior)
end
"""
    make_parameterpoints(pointspecs::NamedTuple)

Create a named tuple of parameter points as name - value pairs based on the given point specifications.

# Arguments
- `pointspecs::NamedTuple`: A named tuple containing the point specifications.

# Returns
A named tuple representing the parameter points, where each field corresponds to a parameter name and its value.

# Examples
```julia
julia> pointspecs = (x = ParameterPointSpec(value = 1), y = ParameterPointSpec(value = 2.5))
julia> make_parameterpoints(pointspecs)
(x = 1, y = 2.5)
"""

function make_parameterpoints(pointspecs::NamedTuple)
    nt = NamedTuple()
    for (k,v) in zip(keys(pointspecs), pointspecs)
        nt = merge(nt, (k => v.value,))
    end
    nt
end