# This file is a part of HS3.jl, licensed under the MIT License (MIT).
"""
    abstract type AbstractHistFactorySpec <: AbstractFunctionalSpec

Abstract type representing the base specification for a HistFactory distribution.
"""
abstract type AbstractHistFactorySpec <: AbstractFunctionalSpec end 

"""
    @with_kw struct HistFactorySpec <: AbstractHistFactorySpec
        axes::NamedTuple{<:Any, <:Tuple{HS3.AxesSpec}}
        samples::NamedTuple
    end

Specification for a HistFactory typed distribution.

# Fields
- `axes::NamedTuple{<:Any, <:Tuple{HS3.AxesSpec}}`: Named tuple containing the axes specifications.
- `samples::NamedTuple`: Named tuple containing the samples.

"""
@with_kw struct HistFactorySpec <: AbstractHistFactorySpec
    axes::NamedTuple{<:Any, <:Tuple{HS3.AxesSpec}}
    samples::NamedTuple
end

#--------------------  HistFactory Sample Specifactions ---------------------------------------------
"""
    @with_kw struct HistFactorySampleSpec <: AbstractHistFactorySpec
        data::AbstractArray{<:Number}
        errors::Union{AbstractArray{<:Number}, Nothing} = nothing 
        modifiers::NamedTuple
    end

Specification for a HistFactory sample.

# Fields
- `data::AbstractArray{<:Number}`: Array containing the sample data.
- `errors::Union{AbstractArray{<:Number}, Nothing}`: Array containing the sample errors, or `nothing` if not applicable.
- `modifiers::NamedTuple`: Named tuple containing any modifiers for the sample.

"""
@with_kw struct HistFactorySampleSpec <: AbstractHistFactorySpec
    data::AbstractArray{<:Number}
    errors::Union{AbstractArray{<:Number}, Nothing} = nothing 
    modifiers::NamedTuple
end

"""
generate_HistFactorySampleSpecs(samples::AbstractArray)

Generate multiple `HistFactorySampleSpec` objects from an array of samples.

# Arguments
- `samples::AbstractArray`: An array of samples.

# Returns
- `sample_specs::NamedTuple`: A named tuple containing the generated `HistFactorySampleSpec` objects.

"""
function generate_HistFactorySampleSpecs(samples)
    sample_specs = NamedTuple()
    for sample in samples
        sample_specs = merge(sample_specs, (_val_content(sample[:name]) => generate_HistFactorySampleSpec(sample),))
    end
    return sample_specs
end

"""
generate_HistFactorySampleSpec(sample::NamedTuple)

Generate a single `HistFactorySampleSpec` object from a named tuple representing a sample.

# Arguments
- `sample::NamedTuple`: A named tuple representing a sample.

# Returns
- `sample_spec::HistFactorySampleSpec`: The generated `HistFactorySampleSpec` object.

"""
function generate_HistFactorySampleSpec(sample)
    modifier_specs = (;)
    for modifier in sample.modifiers
        if haskey(modifier, :parameter) ###for ttW, Parameter sometimes replaces name in that file...
            modifier_specs = merge(modifier_specs, (_val_content(modifier.parameter) => generate_ModifierSpec(modifier),))
        else
            modifier_specs = merge(modifier_specs, (_val_content(modifier.name) => generate_ModifierSpec(modifier),))
        end
    end
    if haskey(sample.data, :errors)
        HistFactorySampleSpec(data = Vector(sample.data.contents), errors = Vector(sample.data.errors), modifiers = modifier_specs)
    else
        HistFactorySampleSpec(data = Vector(sample.data.contents), modifiers = modifier_specs)
    end
end


#-------------------- HistFactory Modifier Specifactions ---------------------------------------------

"""
AbstractModifierSpec <: AbstractHistFactorySpec

An abstract type representing a modifier specification in the HistFactory framework.
This serves as a base type for defining specific modifier specifications.
"""
abstract type AbstractModifierSpec <: AbstractHistFactorySpec end 


"""
ShapesysSpec <: AbstractModifierSpec

A type representing an uncorrelated shape systematics modifier specification in the HistFactory framework.

Fields:
- constraint::Union{Symbol, Nothing}: A constraint associated with the modifier.
- data::NamedTuple: The data associated with the modifier.
"""
@with_kw struct ShapesysSpec <: AbstractModifierSpec 
    constraint::Union{Val{:Gauss}, Val{:Poisson}, Val{:LogNormal}, Val{:Const}, Nothing} = nothing 
    data::Union{NamedTuple, Nothing} = nothing
end

"""
HistosysSpec <: AbstractModifierSpec

A type representing a correlated shape systematics modifier specification in the HistFactory framework.

Fields:
- constraint::Union{Symbol, Nothing}: The constraint associated with the modifier.
- data::NamedTuple: The data associated with the modifier.
"""
@with_kw struct HistosysSpec <: AbstractModifierSpec 
    constraint::Union{Val{:Gauss}, Val{:Poisson}, Val{:LogNormal}, Val{:Const}, Nothing} = nothing 
    data::NamedTuple
end

"""
NormsysSpec <: AbstractModifierSpec

A type representing a normalization systematics modifier specification in the HistFactory framework.

Fields:
- constraint::Union{Symbol, Nothing}: The constraint associated with the modifier.
- data::NamedTuple: The data associated with the modifier.
"""

@with_kw struct NormsysSpec <: AbstractModifierSpec 
    constraint::Union{Val{:Gauss}, Val{:Poisson}, Val{:LogNormal}, Val{:Const}, Nothing} = nothing 
    data::NamedTuple 
end

"""
NormfactorSpec <: AbstractModifierSpec

A type representing a normalization factor modifier specification in the HistFactory framework.
"""
struct NormfactorSpec <: AbstractModifierSpec end

"""
StaterrorSpec <: AbstractModifierSpec

A type representing a statistical error modifier specification in the HistFactory framework.

Fields:
- constraint::Union{Symbol, Nothing}: The constraint associated with the modifier.
"""
@with_kw struct StaterrorSpec <: AbstractModifierSpec 
    constraint::Union{Val{:Gauss}, Val{:Poisson}, Val{:LogNormal}, Val{:Const}, Nothing} = nothing 
end
"""
ShapefactorSpec <: AbstractModifierSpec

A type representing a shape factor modifier specification in the HistFactory framework.

Fields:
- constraint::Union{Symbol, Nothing}: The constraint associated with the modifier.
- data::AbstractArray: The data associated with the modifier.
"""
@with_kw struct ShapefactorSpec <: AbstractModifierSpec 
    constraint::Union{Val{:Gauss}, Val{:Poisson}, Val{:LogNormal}, Val{:Const}, Nothing} = nothing 
    data::AbstractArray = [] 
end

"""
generate_ModifierSpec(nt::NamedTuple) -> AbstractModifierSpec

Generate a modifier specification based on the given NamedTuple.

The type of the modifier specification is determined by the `type` field in the NamedTuple.

Arguments:
- `nt::NamedTuple`: The named tuple containing the modifier specification.

Returns:
- `AbstractModifierSpec`: The generated modifier specification.
"""
function generate_ModifierSpec(nt::NamedTuple)
    tp = nt.type
    nt = NamedTupleTools.delete(nt, (:type, :name, :constraint_name, :parameter)) #for now
    generate_ModifierSpec(tp, nt)
end

"""
generate_ModifierSpec(tp::Type{<:AbstractModifierSpec}, mod::NamedTuple) -> AbstractModifierSpec

Generate a modifier specification of the specified type based on the given NamedTuple.

Arguments:
- `tp::Type{<:AbstractModifierSpec}`: The type of modifier specification to generate.
- `mod::NamedTuple`: The named tuple containing the modifier specification.

Returns:
- `AbstractModifierSpec`: The generated modifier specification.
""" 
generate_ModifierSpec(::Val{:shapefactor}, mod::NamedTuple) = ShapefactorSpec(; mod...)

generate_ModifierSpec(::Val{:shapesys}, mod::NamedTuple) = ShapesysSpec(; mod...)

generate_ModifierSpec(::Val{:histosys}, mod::NamedTuple) = HistosysSpec(; mod...)

generate_ModifierSpec(::Val{:staterror}, mod::NamedTuple) = StaterrorSpec(; mod...)

generate_ModifierSpec(::Val{:normfactor}, mod::NamedTuple) = NormfactorSpec(; mod...)

generate_ModifierSpec(::Val{:normsys}, mod::NamedTuple) = NormsysSpec(; mod...)

