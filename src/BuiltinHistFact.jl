abstract type AbstractHistFactorySpec <: AbstractHS3Spec end 

Parameters.@with_kw struct HistFactorySpec <: AbstractHistFactorySpec
    axes::NamedTuple{<:Any, <:Tuple{HS3.AxesSpec}}
    samples::AbstractArray{AbstractModifierSpec}
end

abstract type AbstractHistFactorySampleSpec end

Parameters.@with_kw struct HistFactorySampleSpec <: AbstractHistFactorySampleSpec
    data::AbstractArray{<:Number}
    errors::Union{AbstractArray{<:Number}, Nothing} = nothing 
    modifiers::NamedTuple
end

abstract type AbstractModifierSpec end 


Parameters.@with_kw struct ShapesysSpec <: AbstractModifierSpec 
    constraint::Union{String, Nothing} = nothing 
    data::AbstractDict{<:Number}
end

Parameters.@with_kw struct HistosysSpec <: AbstractModifierSpec 
    constraint::Union{String, Nothing} = nothing 
    hi::AbstractArray{<:Number}
    lo::AbstractArray{<:Number}
end

Parameters.@with_kw struct NormsysSpec <: AbstractModifierSpec 
    constraint::Union{String, Nothing} = nothing 
    hi::Number
    lo::Number
end

Parameters.@with_kw struct NormfactorSpec <: AbstractModifierSpec end


Parameters.@with_kw struct StaterrorSpec <: AbstractModifierSpec 
    constraint::Union{String, Nothing} = nothing 
    data::AbstractArray = [] 
end

Parameters.@with_kw struct ShapefactorSpec <: AbstractModifierSpec 
    constraint::Union{String, Nothing} = nothing 
    data::AbstractArray = [] 
end

generate_ModifierSpec(::Val{:shapefactor}, mod::NamedTuple) = ShapefactorSpec(; mod...)
generate_ModifierSpec(::Val{:shapesys}, mod::NamedTuple) = ShapesysSpec(; mod...)
generate_ModifierSpec(::Val{:histosys}, mod::NamedTuple) = HistosysSpec(; mod...)
generate_ModifierSpec(::Val{:staterror}, mod::NamedTuple) = StaterrorrSpec(; mod...)
generate_ModifierSpec(::Val{:normfactor}, mod::NamedTuple) = NormfactorSpec(; mod...)
generate_ModifierSpec(::Val{:normsys}, mod::NamedTuple) = NormsysSpec(; mod...)


generate_distspec(::Val{:histfactory_dist}, data::NamedTuple) = HistFactorySpec((samples = generate_HistFactorySampleSpec(data.samples), axes = make_axesspecs(data.axes),))
function generate_distribution(spec::DistributionSpec{:histfactory_dist}) 
    @warn "not yet implemented"
end

function generate_HistFactorySampleSpec(samples::AbstractArray)
    modifiers = NamedTuple()
    for i in samples[:modifiers] 
        modifiers = merge(modifiers, Symbol(i[:name]) => (generate_ModifierSpec(Val(Symbol(i[:type])), filer(x ->x[1] ∉ [:type, :name] ,i))),)
    end
end