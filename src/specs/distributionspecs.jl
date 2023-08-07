# This file is a part of HS3.jl, licensed under the MIT License (MIT).
"""
    abstract type AbstractFunctionalSpec <: AbstractHS3Spec

This abstract type serves as a common interface for various types of functional specifications,
including both distribution and function specifications.
"""
abstract type AbstractFunctionalSpec <: AbstractHS3Spec end


#-------------------- Distribution Specifactions ---------------------------------------------

"""
```julia

    abstract type AbstractDistributionSpec <: AbstractFunctionalSpec

An abstract type representing a distribution specification within the HS3 framework.
"""
abstract type AbstractDistributionSpec <: AbstractFunctionalSpec end


"""
```julia
    struct DistributionSpec{disttype, P<:NamedTuple, T<:NamedTuple} <: AbstractDistributionsSpec

A type representing a distribution specification within the HS3 framework.

# Fields
- `params::P`: A named tuple containing the parameters specific to the distribution.
- `variables::T`: A named tuple containing the variables/variates related to the distribution. 
Depending on the distribution this field is unspecified. In this case the field contains an empty named tuple.

# Type Parameters
- `disttype`: The specific distribution type.
- `P<:NamedTuple`: The type of the `params` field.
- `T<:NamedTuple`: The type of the `variables` field.

# Examples
```julia
spec = DistributionSpec{:gaussian_dist}((mean = 0.0, sigma = 1.0), (x = [1, 2, 3, 4, 5]))
"""


struct DistributionSpec{disttype, P<:NamedTuple} <: AbstractDistributionSpec
    params::P
end

DistributionSpec{disttype}(params::P) where {disttype, P<:NamedTuple} = DistributionSpec{disttype,P}(params)


"""
parse_and_generate_distributionspec(dict::AbstractDict)

Parse a dictionary and delegate the generation of a distribution specification based on the provided data to the corresponding
`generate_distributionspec` function.

Returns a distribution specification of type `DistributionSpec` appropriate to the specific type of distribution. 

# Arguments
- `dict::AbstractDict`: A dictionary containing the necessary parameters for the distribution specification. 
The dictionary should include a `:type` field specifying the distribution type, along with other relevant fields specific to the distribution.

# Examples
```julia
dict = Dict(:type => :gaussian_dist, :mean => 0.0, :sigma => 1.0, :x => [1, 2, 3, 4, 5])
spec = parse_and_generate_distributionspec(dict)
"""

function parse_and_generate_distributionspec(dict::AbstractDict) 
    tp = Val(Symbol(dict[:type]))
    content_dict = filter(entry -> entry[1] != :type, dict)
    generate_distributionspec(tp, _typed_content(content_dict))
end

"""
    generate_distributionspec(tp::Val, data::NamedTuple)

Generate a distribution specification based on the provided type and data.

# Arguments
- `tp::Val`: A `Val` type representing the specific distribution type.
- `data::NamedTuple`: A named tuple containing the necessary data for the distribution specification.

# Returns
A `DistributionSpec` object representing the distribution specification.

# Examples
```julia
data = (mean = 0.0, sigma = 1.0, x = [1, 2, 3, 4, 5])
spec = generate_distributionspec(Val(:gaussian_dist), data)
"""

generate_distributionspec(::Val{:exponential_dist}, data::NamedTuple) = DistributionSpec{:exponential_dist}((c = data.c, x = data.x,))

generate_distributionspec(::Val{:normal_dist}, data::NamedTuple) = generate_distributionspec(Val(:gaussian_dist), data)

generate_distributionspec(::Val{:gaussian_dist}, data::NamedTuple) = DistributionSpec{:gaussian_dist}((mean = data.mean, sigma = data.sigma, x = data.x,))

generate_distributionspec(::Val{:lognormal_dist}, data::NamedTuple) = DistributionSpec{:lognormal_dist}((mean = data.mean, k = data.k, x = data.x,))

generate_distributionspec(::Val{:multinormal_dist}, data::NamedTuple) = DistributionSpec{:multinormal_dist}((covariances = reduce(hcat, data.covariances), mean = data.mean, x = data.x,) )

generate_distributionspec(::Val{:poisson_dist}, data::NamedTuple) = DistributionSpec{:poisson_dist}((mean = data.mean, x = data.x,))

generate_distributionspec(::Val{:argus_dist}, data::NamedTuple) = DistributionSpec{:argus_dist}((resonance = data.resonance, slope = data.slope, power = data.power, x = data.mass,))

generate_distributionspec(::Val{:mixture_dist}, data::NamedTuple) = DistributionSpec{:mixture_dist}((summands = data.summands, coefficients = data.coefficients,)) 

generate_distributionspec(::Val{:product_dist}, data::NamedTuple) = DistributionSpec{:product_dist}((factors = data.factors, )) 

generate_distributionspec(::Val{:uniform_dist}, data::NamedTuple) = DistributionSpec{:uniform_dist}((min = data.min, max = data.max,))

generate_distributionspec(::Val{:polynomial_dist}, data::NamedTuple) = DistributionSpec{:polynomial_dist}((coefficients = data.coefficients, x= data.x, ))

generate_distributionspec(::Val{:CrystalBall_dist}, data::NamedTuple) = DistributionSpec{:CrystalBall_dist}((alpha = data.alpha, x= data.x, m0 = data.m0, n = data.n, sigma = data.sigma, ))

function generate_distributionspec(::Val{:histfactory_dist}, data::NamedTuple)
    HistFactorySpec(; (samples = generate_HistFactorySampleSpecs(data.samples), axes = generate_axes_specs(data.axes),)...)
end


"""
    generate_distributions_specs(dist_array::AbstractArray)
Generate distribution specifications by iterating over an array of distributions and generating individual distribution specs.

# Arguments
- `dist_array`: array containing the distributions.

# Returns
A named tuple containing the generated distribution specifications.

"""
function generate_distributions_specs(dist_array::AbstractArray)
    #dist_specs = NamedTuple{}()
    #for dist in dist_array
    #    dist_name = Symbol(dist.name)
    #    dist_spec = parse_and_generate_distributionspec(dist)
    #    dist_specs = merge(dist_specs, (dist_name => dist_spec,),)
    #end
    #return dist_specs
    names = [Symbol(dist.name) for dist in dist_array]
    specs = [parse_and_generate_distributionspec(dist) for dist in dist_array]
    (; zip(names, specs)...)
end 

#function generate_distributions_specs(dist_array::AbstractArray)
#    chunks = Iterators.partition(dist_array, length((dist_array)) ÷ Threads.nthreads()+1)
#    tasks = map(chunks) do chunk
#        Threads.@spawn _create_distribution_specs(chunk)
#    end
#    chunk_sums = fetch.(tasks)
#    return merge(chunk_sums...)
#end