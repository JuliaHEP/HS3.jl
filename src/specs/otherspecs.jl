# This file is a part of HS3.jl, licensed under the MIT License (MIT).

#-------------------- Uncertainty Specifactions ---------------------------------------------
"""
    @with_kw struct UncertaintySpec{N} <: AbstractUncertaintySpec
Represents a specification of uncertainty of an contents array in the HS3 model as used with binned data sets.

# Fields
- `type::Symbol`: The type of uncertainty. For now only `gaussian_uncertainty` is possible within HS3.
- `sigma::AbstractArray{<:Number, 1}`: An array of standard deviations for each entry of the contents array.
- `correlation::Union{<:Integer, AbstractArray{<:Number, 2}} = 0`: The correlation matrix between 
the entries of the contents array. If no correlation is present, this field should be set to 0 which is also the default behavior. 
If there is correlation, this field should be set to a matrix with the same number of rows and columns as `sigma`.


# Type Parameters
- `N`: Size of the uncertainty array `sigma` and dimensions(N x N) of the `correlation` matrix.


# Examples
```julia
# Create an UncertaintySpec for a variable with Gaussian uncertainty
spec = UncertaintySpec{3}(
    type = :gaussian_uncertainty,
    sigma = [0.1, 0.2, 0.3],
    correlation = [1.0 0.5 0.2;
                   0.5 1.0 0.3;
                   0.2 0.3 1.0]
)

"""
@with_kw struct UncertaintySpec{N} <: AbstractHS3Spec 
    type::Symbol 
    sigma::AbstractArray{<:Number, 1}
    correlation::Union{<:Integer, AbstractArray{<:Number, 2}} = 0
end


"""
    generate_uncertaintyspec(dict::AbstractDict)

Constructs an `UncertaintySpec` object based on the provided abstract dictionary.

# Arguments
- `dict::AbstractDict`: An abstract dictionary containing the specifications for the uncertainty.

# Returns
An `UncertaintySpec` object representing the specified uncertainty.

# Examples
```julia
dict = Dict(:type => "gaussian_uncertainty", :sigma => [0.1, 0.2, 0.3])
spec = generate_uncertaintyspec(dict)
"""
function generate_uncertaintyspec(dict::AbstractDict)
    @assert (dict[:type] == "gaussian_uncertainty") "specified uncertainty type not (yet) supported in HS3"
    N = length(dict[:sigma])
    correlation = haskey(dict, :correlation) ? hcat(dict[:correlation]...) : 0
    nt = NamedTuple{(:type, :sigma, :correlation)}(Symbol(dict[:type]), dict[:sigma], correlation)
    @assert (typeof(size(nt[:correlation])) <: Tuple{M, M} where M<:Integer && (size(nt[:correlation])[1] == N || nt[:correlation] == zeros(1, 1))) "The sizes of the arrays in an uncertainty definition are not correct"
    UncertaintySpec{N}(; nt...)
end


#-------------------- AxesSpecifactions ---------------------------------------------
"""
    @with_kw struct AxesSpec <: AbstractAxesSpec
Struct representing the specification of axes in the HS3 framework.

# Fields
- `max::Union{Number, Nothing} = nothing`: The maximum value of the axis.
- `min::Union{Number, Nothing} = nothing`: The minimum value of the axis.
- `nbins::Union{Integer, Nothing} = nothing`: The number of bins for the axis.

All fields default are optional and default to `nothing`.
"""
@with_kw struct AxesSpec <: AbstractHS3Spec
    max::Union{Number, Nothing} = nothing
    min::Union{Number, Nothing} = nothing 
    nbins::Union{Integer, Nothing} = nothing
end


"""
    generate_axes_specs(axes_array::AbstractArray)

Create axes specifications from an array of axis data.

# Arguments
- `axes_array::AbstractArray`: An array containing axis data.

# Returns
A named tuple representing the axes specifications.

# Example
```julia
axes_data = [...]  # Array of axis data
specs = generate_axes_specs(axes_data)
"""

function generate_axes_specs(axes_array)
    specs = NamedTuple{}()
    for element in axes_array
        temp = NamedTupleTools.delete(_typed_content(element), :name)
        specs = merge(specs, (Symbol(element.name) => AxesSpec(; temp...),),)
    end
    return specs
end


#-------------------- Metadata Specifactions ---------------------------------------------
"""
Metadata contains information about the HS3 file.

# Fields
- `hs3_version::String`: HS3 version number for reference.
- `packages::Vector{Dict{Symbol, String}}`: Optional array of package names and their version numbers used in the creation of this file.
- `authors::Vector{String}`: Optional array of authors, either individual persons or collaborations.
- `publications::Vector{String}`: Optional array of document identifiers of publications associated with this file.
- `description::String`: Optional short abstract/description for this file.
""" 

@with_kw struct MetadataSpec <: AbstractHS3Spec
    hs3_version::String
    packages::AbstractArray = []
    authors::Vector{String} = []
    publications::Vector{String} = []
    description::String = ""
end

"""
generate_metadata(dict)

Generate a MetadataSpec object using the provided dictionary.

# Arguments
- `dict::AbstractDict`: A dictionary containing the metadata information as defined above.

# Returns
A `MetadataSpec` object representing the metadata information.

"""
function generate_metadata_specs(dict::AbstractDict)
    return MetadataSpec(; NamedTuple(dict)...)
end