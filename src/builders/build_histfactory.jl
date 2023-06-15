# This file is a part of HS3.jl, licensed under the MIT License (MIT).
"""
    struct HistfactPDF

The `HistfactPDF` struct represents a probability density function (PDF) defined within the context of a HistFactory model given by LiteHF.jl.

# Fields
- `channel::Tuple`: A tuple representing the channel information associated with the PDF.
- `prior::NamedTuple`: A named tuple representing the prior/modifier information associated with the PDF.

"""
struct HistfactPDF 
    channel::Tuple 
    prior::NamedTuple 
end


"""
    make_dict_from_modifier(modifier) -> Dict

Create a dictionary representation of a given `modifier` object usable for the LiteHF.jl Package.

# Arguments
- `modifier`: The input object.

# Returns
A dictionary representing the input `modifier` object.

"""
make_dict_from_modifier(modifier::NormfactorSpec) = Dict(:type => "normfactor", :data => nothing)

make_dict_from_modifier(modifier::NormsysSpec) = Dict(:type => "normsys", :data => Dict(:hi => modifier.data.hi, :lo => modifier.data.lo))

make_dict_from_modifier(modifier::StaterrorSpec) = Dict(:type => "staterror", :data => modifier.data)

make_dict_from_modifier(modifier::HistosysSpec) = Dict(:type => "histosys", :data => Dict(:hi_data => modifier.data.hi.contents, :lo_data => modifier.data.lo.contents))

make_dict_from_modifier(modifier::ShapesysSpec) = Dict(:type => "shapesys", :data => modifier.data.vals)

make_dict_from_modifier(modifier::ShapefactorSpec) = Dict(:type => "shapefactor", :data => modifier.data)


"""
    make_histfact(histfact_spec::HistFactorySpec, channel_name::Symbol)

Constructs a HistfactPDF object based on the provided HistFactorySpec and channel name using the LiteHF.jl package. 
The specifiactions are first restructured to accomodate for the structure needed by LiteHF.jl which is then used to build a HistFactory object.

# Arguments
- `histfact_spec::HistFactorySpec`: The HistFactorySpec containing the specification for the HistFactory model.
- `channel_name::Symbol`: The name of the channel associated with the HistFactoryPDF.

# Returns
A HistfactPDF object representing the constructed probability density function.

"""
function make_histfact(histfact_spec::HistFactorySpec, channel_name::Symbol)
    sample_names = keys(histfact_spec.samples)
    sample_arr = []
    global_modnames = []
    custom_constraints = []
    for sample in 1:length(histfact_spec.samples)
        temp_dict =  Dict(key => getfield(histfact_spec.samples[sample], key) for key ∈ fieldnames(HistFactorySampleSpec))
        temp_dict[:name] = sample_names[sample]   
        mod_arr = []
        mod_names_in_sample = keys(histfact_spec.samples[sample].modifiers)
        for (mod_name, mod) in zip(keys(histfact_spec.samples[sample].modifiers), histfact_spec.samples[sample].modifiers)
            mod_dict = make_dict_from_modifier(mod)
            mod_dict[:name] = String(mod_name)
            mod_arr = push!(mod_arr, mod_dict)
            #if haskey(mod, :constraint) 
            #    custom_constraints = push!(custom_constraints, (mod_dict[:name], mod.constraint))
            #end
        end
        temp_dict[:modifiers] = mod_arr
        sample_arr  = push!(sample_arr, temp_dict)
        global_modnames = push!(global_modnames, mod_names_in_sample...)
    end
    channel_dict = LiteHF.build_channel(Dict(:samples => sample_arr); misc=Dict())
    channel = LiteHF.build_pyhfchannel((channel_name => channel_dict), global_modnames)
    global_unique = Tuple(unique!(global_modnames))
    #@info channel_name global_unique
    #println([(global_modnames[1][k]) for k in 1:2])
    input_modifiers = [channel[2][k] for k in keys(channel[2])]
    priors = NamedTuple{Tuple(keys(channel[2]))}(LiteHF._prior.(input_modifiers))
    #priors = merge(priors,  NamedTuple{[custom_constraints[2][k] for k in 1:length(custom_constraints)]}(custom_constraints[1]...))
    HistfactPDF(channel, priors)
end

