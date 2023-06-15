# This file is a part of HS3.jl, licensed under the MIT License (MIT).

"""
    ProductDomainSpec <: AbstractHS3Spec

Domain specification representing a product domain in HS3.

# Fields
- `max::Number`: The maximum value of the domain.
- `min::Number`: The minimum value of the domain.
"""
@with_kw struct ProductDomainSpec <: AbstractHS3Spec
    max::Number
    min::Number 
end

"""
    generate_domainspecs(::Val{:product_domain}, arr::AbstractArray)

Generate domain specifications for product domains from an array of dictionaries representing different axes in that domain.
The function iterates over each axis in the array and generates the corresponding domain specification based on the product domain type.

# Arguments
- `arr::AbstractArray`: An array of dictionaries representing the different axes of a domain.

# Returns
A named tuple containing the generated axes specifications of the domain, where the keys are symbols representing the names of the axes.

# Examples
```julia
array = [
    Dict("name" => "domain1", ...),
    Dict("type" => "product_domain", "name" => "domain2", ...)
]

specs = generate_domainspecs(Val(:product_domain), array)
"""

function generate_domainspecs(::Val{:product_domain}, axes::AbstractArray)
    nt = NamedTuple()
    for element in axes
        temp = NamedTuple(filter(entry -> entry[1] != :name, element))
        nt = merge(nt, (Symbol(element.name) => ProductDomainSpec(; temp...),))
    end
    return nt
end

"""
    generate_domains_specs(domain_array::AbstractArray)

Generate domain specifications from an array of domain dictionaries.
The function iterates over each domain dictionary and generates the corresponding domain specification based on its type.
*Note: for now only product domains are part of HS3. Thus, this section might undergo major changes in the future.*

# Arguments
- `domain_array::AbstractArray`: An array of domain dictionaries.

# Returns
A named tuple containing the generated domain specifications, where the keys are symbols representing the names of the domains.

# Examples
```julia
domain_array = [
    Dict("type" => "product_domain", "name" => "domain1", "axes" => [...]),
    Dict("type" => "product_domain", "name" => "domain2", "axes" => [...]),
    ]

specs = generate_domains_specs(domain_array)
"""

function generate_domains_specs(domain_array::AbstractArray)
    specs = NamedTuple()
    for domain in domain_array
        specs = merge(specs, (Symbol(domain.name) => generate_domainspecs(Val(Symbol(domain.type)), domain.axes),))
    end
    specs
end
