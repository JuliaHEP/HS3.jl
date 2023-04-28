# This file is a part of HS3.jl, licensed under the MIT License (MIT).
abstract type AbstractAxesSpec <: AbstractHS3Spec end



Parameters.@with_kw struct AxesSpec <: AbstractAxesSpec 
    max::Union{Number, Nothing} = nothing
    min::Union{Number, Nothing} = nothing 
    nbins::Union{Integer, Nothing} = nothing
end


function make_axesspecs(axes_array::AbstractArray)
    nt = NamedTuple()
    for element in axes_array
        temp = NamedTuple(filter(entry -> entry[1] != :name, element))
        nt = merge(nt, (Symbol(element.name) => AxesSpec(; temp...),),)
    end
    return nt
end
