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
    names = [keys(domain_specs)...]
    values = [LinRange(v.min, v.max, 100) for v in domain_specs]
   (; zip(names, values)...) 
end

"""
    generate_uniform_prior_from_domainspecs(domain_specs::NamedTuple)

Generate a uniform prior distribution from the given domain specifications.
*This function only exist for testing convenience. It will be removed in the future.*

# Arguments
- `domain_specs::NamedTuple`: The domain specifications as a named tuple.

# Returns
A `ValueShapes.NamedTupleDist` object representing the uniform prior distribution.

# Example
```julia
domain_specs = (x = (min = 0, max = 1), y = (min = -1, max = 1))
prior = generate_uniform_prior_from_domainspecs(domain_specs)

"""
function generate_uniform_prior_from_domainspecs(domain_specs::NamedTuple)
    prior = [(Symbol(k), Distributions.Uniform(Float64((v.min)), Float64((v.max)))) for (k, v) in pairs(domain_specs)]
    return ValueShapes.NamedTupleDist(NamedTuple(prior))
end

function generate_uniform_prior_from_domain(domain::NamedTuple)
    tasks = [Threads.@spawn begin
                (Symbol(k), Distributions.Uniform(Float64(first(v)), Float64(last(v))))
            end for (k, v) in pairs(domain)]

    prior_pairs = fetch.(tasks)
    return ValueShapes.NamedTupleDist(NamedTuple(prior_pairs))
end

function generate_ranges_from_domain(domain::NamedTuple)
    prior = [(Symbol(k), ((Float64(first(v)), Float64(last(v))))) for (k, v) in pairs(domain)]
    return NamedTuple(prior)
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
    point_pairs = [(k, Float64(v.value)) for (k, v) in pairs(pointspecs)]
    return NamedTuple(point_pairs)
end


function make_const_parameterpoints(pointspecs::NamedTuple)

    tasks = [Threads.@spawn begin
        (k, Float64(v.value))
    end for (k, v) in pairs(pointspecs) if v.isconst]

    const_pairs = fetch.(tasks) #[(k, Float64(v.value)) for (k, v) in pairs(pointspecs) if v.isconst]
    return NamedTuple(const_pairs)
end

