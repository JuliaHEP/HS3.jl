# This file is a part of HS3.jl, licensed under the MIT License (MIT).
abstract type AbstractUncertaintySpec <: AbstractHS3Spec end

Parameters.@with_kw struct UncertaintySpec{N} <: AbstractUncertaintySpec 
    type::Symbol 
    sigma::AbstractArray{<:Number, 1}
    correlation::Union{<:Integer, AbstractArray{<:Number, 2}} = 0
end

function generate_uncertaintyspec(dict::AbstractDict)
    @assert (dict[:type] == "gaussian_uncertainty") "specified uncertainty type not (yet) supported in HS3"
    nt = NamedTuple()
    N = length(dict[:sigma])
    if haskey(dict, :correlation)
        if dict[:correlation] != 0
            nt = merge(nt, (:correlation => hcat(dict[:correlation]...),))
        else 
            nt = merge(nt, (:correlation => 0,))
        end
        @assert ((typeof(size(nt[:correlation])) <: Tuple{N, N} where N<:Integer && size(nt[:correlation])[1] == N) || nt[:correlation] == 0) "The sizes of the arrays in an uncertainty definition are not correct"
    end
    nt = merge(nt, (:type => Symbol(dict[:type]), :sigma => dict[:sigma],))
    UncertaintySpec{N}(; nt...)
end


"""
Binned Data 
"""

abstract type AbstractDataSpec <: AbstractHS3Spec end

abstract type AbstractBinnedDataSpec <: AbstractDataSpec end 

Parameters.@with_kw struct BinnedDataSpec{N} <: AbstractBinnedDataSpec
    contents::AbstractArray{<:Number, 1} 
    axes::NamedTuple{<:Any, <:Tuple{HS3.AxesSpec}}
    uncertainty::Union{Nothing, UncertaintySpec{N}} = nothing
end

function generate_binneddataspec(dict::AbstractDict)
    content_dict = filter(entry -> entry[1] ∉ [:type, :name], dict)
    content_dict[:axes] = make_axesspecs(dict[:axes])
    if haskey(content_dict, :uncertainty)
        content_dict[:uncertainty] = generate_uncertaintyspec(dict[:uncertainty])
    end
    temp = NamedTuple(content_dict)
    BinnedDataSpec{length(temp.contents)}(; temp...)
end


function make_data(spec::AbstractBinnedDataSpec)
    number_of_bins =  length(spec.contents) 
    min = 0
    max = number_of_bins 
    if spec.axes[1].max !== nothing 
        max = spec.axes[1].max 
    end
    if spec.axes[1].min !== nothing 
        min = spec.axes[1].min 
    end
    step = (max - min)/number_of_bins
    bins = min:step:max
    v = Array(spec.contents)
    StatsBase.Histogram(collect(bins), v)
end



"""
Unbinned Data 
"""

abstract type AbstractUnbinnedDataSpec <: AbstractDataSpec end 

Parameters.@with_kw struct UnbinnedDataSpec{N} <: AbstractUnbinnedDataSpec 
    entries::AbstractArray{<:Number} 
    entries_errors::Union{Nothing, AbstractArray{<:Number}} = nothing 
    axes::NamedTuple{<:Any, <:NTuple{N, AxesSpec}}
    weights::Union{AbstractArray{<:Number}, Nothing}  = nothing
end

function generate_unbinneddataspec(dict::AbstractDict)
    content_dict = filter(entry -> entry[1] ∉ [:type, :name], dict)
    content_dict[:axes] = make_axesspecs(dict[:axes])
    if haskey(dict, :entries_errors)
        @assert size(dict[:entries_errors]) == size(dict[:entries]) "Size/shape of entries_errors in unbinned data not correct"
        content_dict[:entries_errors] = reduce(vcat, dict[:entries_errors]')
    end 
    if haskey(dict, :weights)
        @assert length(dict[:weights]) == length(dict[:entries]) "Size of weights in unbinned data not correct"
        content_dict[:weights] = Array(dict[:weights])
    end
    content_dict[:entries] = reduce(vcat, (dict[:entries]))
    print(typeof(content_dict[:entries]))
    temp = NamedTuple(content_dict)
    UnbinnedDataSpec{length(temp.entries)}(; temp...)
end

abstract type AbstractPointDataSpec <: AbstractDataSpec end

Parameters.@with_kw struct PointDataSpec <: AbstractPointDataSpec 
    value::Number
    error::Union{Number, Nothing} = nothing 
end

function generate_pointdataspec(dict::AbstractDict)
    content_dict = filter(entry -> entry[1] ∉ [:type, :name], dict)
    PointDataSpec(; NamedTuple(content_dict)...)
end

function make_dataspecs(dict::AbstractArray)
    nt = NamedTuple()
    for data_set in dict 
        if data_set.type == "unbinned"
            nt = merge(nt, (Symbol(data_set.name) => generate_unbinneddataspec(data_set),))
        elseif data_set.type == "binned"
            nt = merge(nt, (Symbol(data_set.name) => generate_binneddataspec(data_set),)) 
        elseif data_set.type == "point"
            nt = merge(nt, (Symbol(data_set.name) => generate_pointdataspec(data_set),))
        else
            @warn "Specified type of data $(data_set.type) not part of HS3"
        end
    end
    nt
end

#BinnedDataSpec{N::Integer}(contents, axes::P) where {P<:AbstractDict} = BinnedDataSpec{N}(contents, make_axesspecs(axes))
#function generate_dataspec(::Val{:unbinned}, data::NamedTuple) 
#    temp_weights = ones(Float64, length(data.entries)) 
#    temp_axes = make_axesspecs(collect(data.axes))
#    temp_error = []
#    if haskey(data, :weights) temp_weights = collect(data.weights) end 
#    if haskey(data, :entries_errors) 
#        temp_error = collect(data.entries_errors)
#        @assert size(collect(temp_error)) == size(collect(data.entries))
#    end 
#    @assert length(temp_weights) == length(collect(data.entries)) 
#    BuiltinDataSpec{:unbinned}((entries = collect(data.entries), axes = temp_axes, weights = temp_weights, entries_error = temp_error,))
#end
#
#function generate_dataspec(::Val{:binned}, data::NamedTuple) 
#    temp_axes = make_axesspecs(collect(data.axes))
#    temp_uncert = nothing 
#    if haskey(data, :uncertainty)
#        temp_uncert = (type = data.uncertainty.type, sigma = data.uncertainty.sigma,)
#        print(temp_uncert.type)
#        @assert temp_uncert.type == Val(Symbol("gaussian_uncertainty")) "For now only Gaussian uncertainties are supoorted for binned data in HS3"
#        @assert length(temp_uncert.sigma) == length(data.contents) "In binned data the sigma of uncertainty needs to be of same length as the contents array"
#        if haskey(data.uncertainty, :correlation)
#            temp_uncert = merge(temp_uncert, (correlation = data.uncertainty.correlation,))
#        else
#            temp_uncert = merge(temp_uncert, (correlation = 0,))
#        end
#    end
#    return BuiltinDataSpec{:binned}((axes = temp_axes, contents = data.contents, uncertainty = temp_uncert, ))
#end
#
#function generate_dataspec(::Val{:point}, data::NamedTuple)
#    haskey(data, :error) ? temp_error = data.error : temp_error = 0
#    @assert typeof(data.value) <: Number
#    BuiltinDataSpec{:point}((value = data.value, error = temp_error,))
#end
#
#
#function make_dataspecs(data_array::AbstractArray)
#    specs = NamedTuple()
#    for data in 1:length(data_array)
#        specs = merge(specs, (Symbol(data_array[data].name) => AbstractDataSpec(data_array[data]),))
#    end
#    return specs
#end
#
#export make_dataspecs
#
#
#
#
#function generate_hist(spec::BuiltinDataSpec{:binned})
#    step = (spec.params.axes[1].max - spec.params.axes[1].min)/spec.params.axes[1].nbins
#    bins = spec.params.axes[1].min:step:spec.params.axes[1].max
#    #print(typeof(bins))
#    #println(spec.params.axes[1].max)
#    #println(collect(bins))
#    v = Array(spec.params.contents)
#    #println(v)
#    #StatsBase.fit(StatsBase.Histogram, spec.params.contents, bins)
#    StatsBase.Histogram(collect(bins), v)
#    #return hist
#    #StatsBase.Histogram([spec.params.coordinates..., spec.params.observables[1].max], [spec.params.counts...])
#end
#
#export generate_hist
#
#function make_hists(specs::NamedTuple)
#    hists = []
#    for spec in specs 
#        hists = push!(hists, generate_hist(spec))
#    end
#    return (; zip(keys(specs), hists)...)
#end
#
#export make_hists
#