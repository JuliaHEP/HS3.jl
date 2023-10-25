# This file is a part of HS3.jl, licensed under the MIT License (MIT).
"""
    abstract type AbstractDataSpec <: AbstractHS3Spec

This abstract type serves as a common interface for various types of data specifications within HS3..
"""
abstract type AbstractDataSpec <: AbstractHS3Spec end

#-------------------- Binned Data Specifactions ---------------------------------------------
"""
    @with_kw struct BinnedDataSpec{N} <: AbstractBinnedDataSpec
        contents::AbstractArray{<:Number, 1}
        axes::NamedTuple{<:Any, <:Tuple{HS3.AxesSpec}}
        uncertainty::Union{Nothing, UncertaintySpec{N}} = nothing
    end

Specification for binned data within HS3.

# Fields
- `contents::AbstractArray{<:Number, 1}`: Array of binned data values.
- `axes::NamedTuple{<:Any, <:Tuple{HS3.AxesSpec}}`: NamedTuple of axes specifications.
- `uncertainty::Union{Nothing, UncertaintySpec{N}}`: Uncertainty specification, defaulting to `nothing`.

# Type Parameters
- `N`: Size of the contents array.
"""
@with_kw struct BinnedDataSpec{N} <: AbstractDataSpec
    contents::AbstractArray{<:Number, 1} 
    axes::NamedTuple#NamedTuple{<:Any, <:Tuple{HS3.AxesSpec}}
    uncertainty::Union{Nothing, UncertaintySpec{N}} = nothing
end


"""
    generate_binneddataspec(dict::AbstractDict)

Constructs a `BinnedDataSpec` object based on the provided abstract dictionary.

# Arguments
- `dict::AbstractDict`: An abstract dictionary containing the specifications for the binned data.

# Returns
A `BinnedDataSpec` object representing the specified binned data.

"""
function generate_dataspec(dict::AbstractDict, ::Val{:binned})
    content_dict = Dict(filter(entry -> entry[1] ∉ [:type, :name], dict))
    content_dict[:axes] = generate_axes_specs(dict[:axes])
    haskey(content_dict, :uncertainty) ? content_dict[:uncertainty] = generate_uncertaintyspec(dict[:uncertainty]) : nothing
    temp = NamedTuple(content_dict)
    BinnedDataSpec{length(content_dict[:contents])}(; temp...)
end


#-------------------- Unbinned Data Specifactions ---------------------------------------------
"""
    struct UnbinnedDataSpec{N} <: AbstractDataSpec

Specification for unbinned data within the HS3 framework.

# Fields
- `entries::AbstractArray{<:Number}`: Array of data entries.
- `entries_errors::Union{Nothing, AbstractArray{<:Number}}`: Array of error values for the data entries (optional).
- `axes::NamedTuple{<:Any, <:NTuple{N, AxesSpec}}`: NamedTuple representing the axes specifications.
- `weights::Union{AbstractArray{<:Number}, Nothing}`: Array of weights for the data entries (optional).

# Type Parameters
- `N`: Size of the contents array.
"""
@with_kw struct UnbinnedDataSpec <: AbstractDataSpec 
    entries::AbstractArray{<:Number} 
    entries_errors::Union{Nothing, AbstractArray{<:Number}} = nothing 
    axes::NamedTuple#NamedTuple{<:Any, <:NTuple{N, AxesSpec}}
    weights::Union{AbstractArray{<:Number}, Nothing}  = nothing
end

"""
    generate_unbinneddataspec(dict::AbstractDict)

Constructs an `UnbinnedDataSpec` object based on the provided abstract dictionary.

# Arguments
- `dict::AbstractDict`: An abstract dictionary containing the specifications for the unbinned data.

# Returns
An `UnbinnedDataSpec` object representing the specified unbinned data.

# Examples
```julia
dict = Dict(
    :type => "unbinned_data",
    :entries => [[1.0], [2.0], [3.0]],
    :entries_errors => [[0.1], [0.2], [0.3]],
    :axes => NamedTuple{Tuple{AxesSpec}}((:x => AxesSpec(:name => "x"))), TODO: check
    :weights => [0.5, 0.3, 0.2]
)

spec = generate_unbinneddataspec(dict)
"""
function generate_dataspec(dict::AbstractDict, ::Val{:unbinned})
    content_dict = filter(entry -> entry[1] ∉ [:type, :name], dict)
    content_dict[:axes] = generate_axes_specs(dict[:axes])
    if haskey(dict, :entries_errors)
        @assert size(dict[:entries_errors]) == size(dict[:entries]) "Size/shape of entries_errors in unbinned data not correct"
        content_dict[:entries_errors] = reduce(vcat, dict[:entries_errors]')
    end 
    if haskey(dict, :weights)
        @assert length(dict[:weights]) == length(dict[:entries]) "Size of weights in unbinned data not correct"
        content_dict[:weights] = Array(dict[:weights])
    end
    content_dict[:entries] = reduce(vcat, (dict[:entries]))
    temp = NamedTuple(content_dict)
    UnbinnedDataSpec(; temp...)
end

#-------------------- Point Data Specifactions ---------------------------------------------
"""
    struct PointDataSpec <: AbstractDataSpec

Specification for point data within the HS3 framework.

# Fields
- `value::Number`: Value of the data point.
- `error::Union{Number, Nothing}`: Error value for the data point (optional).

"""
@with_kw struct PointDataSpec <: AbstractDataSpec 
    value::Number
    error::Union{Number, Nothing} = nothing 
end

"""
    function generate_pointdataspec(dict::AbstractDict) -> PointDataSpec

Generate a `PointDataSpec` instance based on the provided dictionary.

# Arguments
- `dict::AbstractDict`: Dictionary containing the specifications for the point data.

# Returns
An `PointDataSpec` object representing the specified unbinned data.

"""
function generate_dataspec(dict::AbstractDict, ::Val{:point})
    content_dict = filter(entry -> entry[1] ∉ [:type, :name], dict)
    PointDataSpec(; NamedTuple(content_dict)...)
end


#-------------------- Automation ---------------------------------------------

"""
    generate_data_specs(array::AbstractArray)

Generate data specifications from an array of dictionaries representing different data sets.
The function iterates over each data set in the array and generates the corresponding data specification based on its type.

# Arguments
- `array::AbstractArray`: An array of dictionaries representing different data sets.

# Returns
A named tuple containing the generated data specifications, where the keys are symbols representing the names of the data sets.

# Examples
```julia
array = [
    Dict("type" => "unbinned", "name" => "dataset1", ...),
    Dict("type" => "binned", "name" => "dataset2", ...),
    Dict("type" => "point", "name" => "dataset3", ...)
]

specs = generate_data_specs(array)

"""
function generate_data_specs(array::AbstractArray)
    #data_specs = NamedTuple()
    names = [Symbol(data_set.name) for data_set in array]
    specs = [generate_dataspec(data_set, Val(Symbol(data_set.type))) for data_set in array]
    #for data_set in array 

        #if data_set.type == "unbinned"
        #    data_specs = merge(data_specs, (Symbol(data_set.name) => generate_unbinneddataspec(data_set),))
        #elseif data_set.type == "binned"
        #    data_specs = merge(data_specs, (Symbol(data_set.name) => generate_binneddataspec(data_set),)) 
        #elseif data_set.type == "point"
        #    data_specs = merge(data_specs, (Symbol(data_set.name) => generate_pointdataspec(data_set),))
        #else
        #    @warn "Specified type of data $(data_set.type) not part of HS3"
        #end
    #end
    return (; zip(names, specs)...)
end