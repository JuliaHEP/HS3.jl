abstract type AbstractHistFactorySpec <: AbstractFunctionalSpec end 

abstract type AbstractModifierSpec <: AbstractHistFactorySpec end 


Parameters.@with_kw struct HistFactorySpec <: AbstractHistFactorySpec
    axes::NamedTuple{<:Any, <:Tuple{HS3.AxesSpec}}
    samples::NamedTuple
end

abstract type AbstractHistFactorySampleSpec end

Parameters.@with_kw struct HistFactorySampleSpec <: AbstractHistFactorySampleSpec
    data::AbstractArray{<:Number}
    errors::Union{AbstractArray{<:Number}, Nothing} = nothing 
    modifiers::NamedTuple
end

##########

Parameters.@with_kw struct ShapesysSpec <: AbstractModifierSpec 
    constraint::Union{String, Nothing} = nothing 
    data::NamedTuple
end

Parameters.@with_kw struct HistosysSpec <: AbstractModifierSpec 
    constraint::Union{String, Nothing} = nothing 
    data::NamedTuple
end

Parameters.@with_kw struct NormsysSpec <: AbstractModifierSpec 
    constraint::Union{String, Nothing} = nothing 
    data::NamedTuple #{Symbol, <:Number}
end

struct NormfactorSpec <: AbstractModifierSpec end


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

function generate_ModifierSpec(nt::NamedTuple)
    tp = nt.type
    #print(nt)
    nt = NamedTupleTools.delete(nt, (:type, :name))
    generate_ModifierSpec(tp, nt)
end


function generate_distspec(::Val{:histfactory_dist}, data::NamedTuple)
    HistFactorySpec(; (samples = generate_HistFactorySampleSpec(data.samples), axes = make_axesspecs(data.axes),)...)
end
function generate_distribution(spec::DistributionSpec{:histfactory_dist}) 
    @warn "not yet implemented"
end

function generate_HistFactorySampleSpec(samples::AbstractArray)
    sample_specs = NamedTuple()
    for sample in samples
        sample_specs = merge(sample_specs, (Symbol(_val_content(sample.name)) => generate_HistFactorySampleSpec(sample),))
    end
    sample_specs
end

function generate_HistFactorySampleSpec(sample::NamedTuple)
    mod_specs = (;)
    println(sample[:modifiers])
    for mod in sample[:modifiers]
        mod_specs = merge(mod_specs, (Symbol(_val_content(mod.name)) => generate_ModifierSpec(mod),))
    end
    HistFactorySampleSpec(; (data = Vector(sample.data.contents), modifiers = mod_specs)...)
end

make_dict_from_mod(mod::NormfactorSpec) = Dict(:type => "normfactor", :data => nothing)

make_dict_from_mod(mod::NormsysSpec) = Dict(:type => "normsys", :data => Dict(:hi => mod.data.hi, :lo => mod.data.lo))

make_dict_from_mod(mod::StaterrorSpec) = Dict(:type => "staterror", :data => mod.data)

make_dict_from_mod(mod::HistosysSpec) = Dict(:type => "histosys", :data => Dict(:hi_data => mod.data.hi.contents, :lo_data => mod.data.lo.contents))

make_dict_from_mod(mod::ShapesysSpec) = Dict(:type => "shapesys", :data => mod.data.vals)

make_dict_from_mod(mod::ShapefactorSpec) = Dict(:type => "shapefactor", :data => mod.data)


struct HistfactPDF 
    channel::Tuple 
    prior::NamedTuple 
end

function make_histfact(histfact_spec::HistFactorySpec, channel_name::Symbol)
    #hist_fact_dict = Dict(key=>getfield(histfact_spec, key) for key ∈ fieldnames(HistFactorySpec))
    sample_names = keys(histfact_spec.samples)
    sample_arr = []
    global_modnames = []
    for sample in 1:length(histfact_spec.samples)
        temp_dict =  Dict(key => getfield(histfact_spec.samples[sample], key) for key ∈ fieldnames(HistFactorySampleSpec))
        temp_dict[:name] = sample_names[sample]   
        mod_arr = []
        mod_names_in_sample = keys(histfact_spec.samples[sample].modifiers)
        for (mod_name, mod) in zip(keys(histfact_spec.samples[sample].modifiers), histfact_spec.samples[sample].modifiers)
            mod_dict = make_dict_from_mod(mod)
            mod_dict[:name] = String(mod_name)
            mod_arr = push!(mod_arr, mod_dict)
        end
        temp_dict[:modifiers] = mod_arr
        sample_arr  = push!(sample_arr, temp_dict)
        global_modnames = push!(global_modnames, mod_names_in_sample...)
       #Dict(hist_fact_dict[:samples][1] => getfield(histfact_spec, key) for key ∈ fieldnames(HistFactorySampleSpec))
    end
    #global_modnames
    channel_dict = LiteHF.build_channel(Dict(:samples => sample_arr); misc=Dict())
    channel = LiteHF.build_pyhfchannel((channel_name => channel_dict), global_modnames)
    global_unique = Tuple(unique!(global_modnames))
    input_modifiers = [channel[2][k] for k in global_unique]
    priors = NamedTuple{global_unique}(LiteHF._prior.(input_modifiers))
    HistfactPDF(channel, priors)
end

