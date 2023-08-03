# This file is a part of HS3.jl, licensed under the MIT License (MIT).
"""
    struct LikelihoodSpec <: AbstractHS3Spec

Specification for likelihood-related data.

This struct represents the specification for likelihoods, including the observed data, the distributions used for modeling, and optional auxiliary distributions.

Fields:
- `data`: An abstract array of strings representing the observed data.
- `distributions`: An abstract array of strings representing the distributions used for modeling.
- `aux_distributions`: An optional abstract array of strings representing auxiliary distributions. If not provided, it defaults to `nothing`.

Example:
```julia
spec = LikelihoodSpec(
    data = ["obs1", "obs2"],
    distributions = ["dist1", "dist2"],
    aux_distributions = ["aux_dist1", "aux_dist2"]
)

"""
@with_kw struct LikelihoodSpec <: AbstractHS3Spec
    data::Union{AbstractArray{String}, Nothing} = nothing
    distributions::Union{AbstractArray{String}, Nothing} = nothing
    aux_distributions::Union{AbstractArray{String}, Nothing} = nothing
end

"""
    generate_likelihoodspec(dict::AbstractDict)

Generate a `LikelihoodSpec` from a dictionary.

This function takes a dictionary `dict` representing the specification for a likelihood and generates a `LikelihoodSpec` object.
The dictionary should contain the following keys:
- `:data`: An array-like object containing the observed data as strings.
- `:distributions`: An array-like object containing the distributions used for modeling as strings.
- `:name` (optional): A string representing the name of the likelihood (ignored by this function).

The function filters out the `:name` key from the dictionary and checks that the lengths of `:data` and `:distributions` are the same. If the lengths are not equal, an `@assert` error is raised.

Example:
```julia
dict = Dict(
    :name => "my_likelihood",
    :data => ["obs1", "obs2"],
    :distributions => ["dist1", "dist2"]
)
spec = generate_likelihoodspec(dict)
"""
function generate_likelihoodspec(dict::AbstractDict)
    content_dict = filter(entry -> entry[1] != :name, dict)
    #@assert length(content_dict[:data]) == length(content_dict[:distributions]) "Likelihood $(dict[:name]) not defined correctly."
    LikelihoodSpec(; NamedTuple(content_dict)...)
end


"""
    generate_likelihoods_specs(arr::AbstractArray)

Generate multiple likelihood specifications from an array.

This function takes an array `arr` containing elements representing multiple likelihood specifications. It iterates over each element in the array, generates a `LikelihoodSpec` for each element, and stores them in a named tuple. The `name` field in each element is used as the key in the named tuple.

# Arguments
- `arr::AbstractArray`: An array containing elements representing multiple likelihood specifications.

# Returns
- `NamedTuple`: A named tuple containing the generated likelihood specifications.

# Example
```julia
dict1 = Dict(
    :name => "likelihood1",
    :data => ["obs1", "obs2"],
    :distributions => ["dist1", "dist2"]
)
dict2 = Dict(
    :name => "likelihood2",
    :data => ["obs3", "obs4"],
    :distributions => ["dist3", "dist4"]
)
arr = [dict1, dict2]
specs = generate_likelihoods_specs(arr)

"""
function generate_likelihoods_specs(arr::AbstractArray)
    names = [Symbol(element.name) for element in arr]
    specs = [generate_likelihoodspec(element) for element in arr]
    #for element in arr
    #    specs = merge(specs, (Symbol(element.name) => generate_likelihoodspec(element),))
    #end
    #specs
    (; zip(names, specs)...)
end
