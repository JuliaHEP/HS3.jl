# This file is a part of HS3.jl, licensed under the MIT License (MIT).

"""
    ProductDomainSpec <: AbstractHS3Spec

Domain specification representing a product domain in HS3.

# Fields
- `max::Number`: The maximum value of the domain.
- `min::Number`: The minimum value of the domain.
"""
@with_kw struct ProductDomainSpec{T <: Real, S <: Real} <: AbstractHS3Spec
    max::T
    min::S 
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
    #nt = (;)
    names = [Symbol(element.name) for element in axes]
    specs = [ProductDomainSpec(element.max, element.min) for element in axes]
    #for element in axes
    #    #if !occursin("gamma_shape_stat", element.name)
    #     #   if occursin("gamma_stat_CombReco_bin", element.name)
    #            bin_number = match(r"bin_(\d+)$", element.name).captures[1]
    #            new_name = "staterror_bin$(bin_number)"
    #            temp = NamedTuple(filter(entry -> entry[1] != :name, element))
    #            nt = merge(nt, (Symbol(new_name) => ProductDomainSpec(temp...),))
    #        else
    #            temp = NamedTuple(filter(entry -> entry[1] != :name, element))
    #            nt = merge(nt, (Symbol(element.name) => ProductDomainSpec(temp...),))
    #        end
    #    end
    #end
    return (; zip(names, specs)...)
end

#function generate_domainspecs(::Val{:product_domain}, axes::AbstractArray)
#    chunks = Iterators.partition(axes, length((axes)) ÷ Threads.nthreads()+1)
#    tasks = map(chunks) do chunk
#        Threads.@spawn _create_domain_axes(chunk)
#    end
#    chunk_sums = fetch.(tasks)
#    return merge(chunk_sums...)
#end

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
    names = [Symbol(domain.name) for domain in domain_array]
    specs = [generate_domainspecs(Val(Symbol(domain.type)), domain.axes) for domain in domain_array]
    #for domain in domain_array
    #    specs = merge(specs, (Symbol(domain.name) => generate_domainspecs(Val(Symbol(domain.type)), domain.axes),))
    #end
    #println("here")
    return (; zip(names, specs)...)
end
