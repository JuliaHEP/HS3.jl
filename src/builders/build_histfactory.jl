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
make_dict_from_modifier(modifier::NormfactorSpec, sample::HistFactorySampleSpec) = Dict(:type => "normfactor", :data => nothing)

make_dict_from_modifier(modifier::NormsysSpec, sample::HistFactorySampleSpec) = Dict(:type => "normsys", :data => Dict(:hi => modifier.data.hi, :lo => modifier.data.lo))

function make_dict_from_modifier(modifier::StaterrorSpec, sample::HistFactorySampleSpec) 
    @assert sample.errors !== nothing
    Dict(:type => "staterror", :data => sample.errors/sum(sample.data))
end
make_dict_from_modifier(modifier::HistosysSpec, sample::HistFactorySampleSpec) = Dict(:type => "histosys", :data => Dict(:hi_data => modifier.data.hi.contents, :lo_data => modifier.data.lo.contents))

function make_dict_from_modifier(modifier::ShapesysSpec, sample::HistFactorySampleSpec)
    if modifier.data === nothing 
        @assert modifier.constraint == Val{:Const}()
        data = zeros(length(sample.data))
    else 
        data = modifier.data.vals 
    end
    Dict(:type => "shapesys", :data => data)
end 

make_dict_from_modifier(modifier::ShapefactorSpec, sample::HistFactorySampleSpec) = Dict(:type => "shapefactor", :data => modifier.data)
#TODO check this

"""
make_modifier_dict(sample::HistFactorySampleSpec)

Create a modifier dictionary for LiteHF.jl and constraint dictionary from a HistFactorySampleSpec object.

# Arguments
- `sample`: A HistFactorySampleSpec object containing modifiers.

# Returns
- `mod_arr`: An array of modifier dictionaries.
- `constraints`: A dictionary mapping modifier names to their constraints, if specified.

"""
function make_modifier_dict(sample::HistFactorySampleSpec)
    constraints = NamedTuple()
    mod_arr = []
    for (mod_name, mod) in zip(keys(sample.modifiers), sample.modifiers)
        mod_dict = make_dict_from_modifier(mod, sample)
        mod_dict[:name] = String(mod_name)
        mod_arr = push!(mod_arr, mod_dict)
        if :constraint ∈ fieldnames(typeof(mod)) && mod.constraint !== nothing
            constraints = merge(constraints, [Symbol(mod_dict[:name]) => mod.constraint,])
        end 
    end
    mod_arr, constraints
end


"""
_make_constraint(::Val{distribution_type})

Create and return a constraint object based on the specified distribution type.

# Arguments
- `distribution_type`: A `Val` type specifying the distribution type.

# Returns
- A constraint object corresponding to the specified distribution type.

# Supported Distribution Types
- `:Gauss`: Gaussian distribution (Normal distribution)
- `:Poisson`: Poisson distribution
- `:Const`: Constant distribution
- `:LogNormal`: Log-normal distribution

"""
### TODO: This needs improvement
_make_constraint(::Val{:Gauss}) = Distributions.Normal()

_make_constraint(::Val{:Poisson}) = LiteHF.RelaxedPoisson(1.)

_make_constraint(::Val{:Const}) = LiteHF.FlatPrior(0, 5)

_make_constraint(::Val{:LogNormal}) = Distributions.LogNormal()


"""
make_sample_dict(samples)

Create a sample dictionary and constraint dictionary from a collection of samples.

# Arguments
- `samples`: A collection of samples.

# Returns
- `sample_arr`: An array of sample dictionaries.
- `constraints`: A dictionary mapping modifier names to their constraints, if specified.

"""
function make_sample_dict(samples)
    constraints = NamedTuple()
    sample_arr = []
    for (sample_name, sample) in zip(keys(samples), samples)
        sample_dict =  Dict(key => getfield(sample, key) for key ∈ fieldnames(HistFactorySampleSpec))
        sample_dict[:name] = sample_name   
        mod_arr, temp_constraints = make_modifier_dict(sample)
        constraints = merge(constraints, temp_constraints)
        sample_dict[:modifiers] = mod_arr
        sample_arr  = push!(sample_arr, sample_dict)
    end
    sample_arr, constraints
end


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
    sample_arr, constraints = make_sample_dict(histfact_spec.samples)
    channel_dict = LiteHF.build_channel(Dict(:samples => sample_arr); misc=Dict())
    modifier_names = unique!(reduce(vcat, [channel_dict[i].modifier_names for i in keys(channel_dict)]))
    channel = LiteHF.build_pyhfchannel((channel_name => channel_dict), modifier_names)
    @assert length(modifier_names) == length(keys(channel[2]))


    constraint_terms = NamedTuple()
    for name in unique!(modifier_names)
        if name in keys(constraints)
            constraint_terms = merge(constraint_terms, [name => _make_constraint(constraints[name])])
        else
            constraint_terms = merge(constraint_terms, [name => LiteHF._prior(channel[2][name])])
        end
    end
    HistfactPDF(channel, constraint_terms)
end

