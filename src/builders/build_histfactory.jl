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
    order::Vector
    customs::NamedTuple
end


"""
    make_dict_from_modifier(modifier) -> Dict

Create a dictionary representation of a given `modifier` object usable for the LiteHF.jl Package.

# Arguments
- `modifier`: The input object.

# Returns
A dictionary representing the input `modifier` object.

"""
make_dict_from_modifier(modifier::NormfactorSpec, sample::HistFactorySampleSpec, abs_staterror_data::Number) = Dict(:type => "normfactor", :data => nothing)

make_dict_from_modifier(modifier::NormsysSpec, sample::HistFactorySampleSpec, abs_staterror_data::Number) = Dict(:type => "normsys", :data => Dict(:hi => modifier.data.hi, :lo => modifier.data.lo))

function make_dict_from_modifier(modifier::StaterrorSpec, sample::HistFactorySampleSpec, abs_staterror_data::Number) 
    @assert sample.errors !== nothing
    Dict(:type => "staterror", :data =>  sample.errors * abs_staterror_data)
end
make_dict_from_modifier(modifier::HistosysSpec, sample::HistFactorySampleSpec, abs_staterror_data::Number) = Dict(:type => "histosys", :data => Dict(:hi_data => modifier.data.hi.contents, :lo_data => modifier.data.lo.contents))

function make_dict_from_modifier(modifier::ShapesysSpec, sample::HistFactorySampleSpec, abs_staterror_data::Number)
    if modifier.data === nothing 
        data = zeros(length(sample.data))
    else 
        data = modifier.data.vals 
    end
    Dict(:type => "shapesys", :data => data)
end 

make_dict_from_modifier(modifier::ShapefactorSpec, sample::HistFactorySampleSpec, abs_staterror_data::Number) = Dict(:type => "shapefactor", :data => modifier.data)

make_dict_from_modifier(modifier::CustomModifierSpec, sample::HistFactorySampleSpec, abs_staterror_data::Number) = Dict(:type => "custom", :data => nothing)

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
function make_modifier_dict(sample::HistFactorySampleSpec, abs_staterror_data::Number, function_specs::NamedTuple)
    constraints = NamedTuple()
    mod_arr = []
    custom = (;)
    for (mod_name, mod) in zip(keys(sample.modifiers), sample.modifiers)
        if typeof(mod) != CustomModifierSpec 
            mod_dict = make_dict_from_modifier(mod, sample, abs_staterror_data)
            mod_dict[:name] = String(mod_name)
            mod_arr = push!(mod_arr, mod_dict)
        else
            sorted = topological_sort(function_specs)
            func = make_functional(function_specs[Symbol(mod_name)], sorted)
            #TODO: Account for non generic functions
            parameter_names = _val_content(collect(function_specs[Symbol(mod_name)].params.var))
            for i in function_specs 
                if typeof(i) == FunctionSpec{:generic_function}
                    parameter_names = append!(parameter_names, _val_content(collect(i.params.var)))
                end
            end
            #@info parameter_names
            mod_dict = Dict(:type=> "custom", :data =>(params -> func(params), parameter_names) , :name => mod_name)
            mod_arr = push!(mod_arr, mod_dict)
            custom = merge(custom, (mod_name => parameter_names,))
        end
        if :constraint ∈ fieldnames(typeof(mod)) && mod.constraint !== nothing
            constraints = merge(constraints, [Symbol(mod_dict[:name]) => mod.constraint,])
        end 
    end
    mod_arr, constraints, custom
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

_make_constraint(::Val{:Gauss}, modifier::LiteHF.Normsys) =  Distributions.Normal(0, 1.)

_make_constraint(::Val{:Gauss}, modifier::LiteHF.Histosys) = Distributions.Normal(0, 1.)

_make_constraint(::Val{:Gauss}, modifier::LiteHF.Shapesys) = Distributions.Normal(0, 1)

_make_constraint(::Val{:Gauss}, modifier::LiteHF.Staterror) = Distributions.Normal(0, modifier.σ)

_make_constraint(::Val{:Poisson}, modifier::LiteHF.Normsys) = LiteHF.RelaxedPoisson(1.)

_make_constraint(::Val{:Poisson}, modifier::LiteHF.Histosys) = LiteHF.RelaxedPoisson(1.)

_make_constraint(::Val{:Poisson}, modifier::LiteHF.Shapesys) = LiteHF.RelaxedPoisson(modifier.moddata^2)

_make_constraint(::Val{:Poisson}, modifier::LiteHF.Staterror) = LiteHF.RelaxedPoisson(modifier.σ^2)

_make_constraint(::Val{:Const}, modifier::LiteHF.AbstractModifier) = LiteHF.FlatPrior(1, 1) # TODO: this needs improvement

_make_constraint(::Any, ::Any) = @error "Not a valid entry and thus not defined"

#_make_constraint(::Val{:LogNormal} ) = Distributions.LogNormal()


function _calc_staterror_data(samples::NamedTuple)
    count = 0.0 
    for sample in samples
        if any(mod -> typeof(mod) == StaterrorSpec, sample.modifiers)
            count += sum(sample.data)
        end
    end
    count
end
"""
make_sample_dict(samples)

Create a sample dictionary and constraint dictionary from a collection of samples.

# Arguments
- `samples`: A collection of samples.

# Returns
- `sample_arr`: An array of sample dictionaries.
- `constraints`: A dictionary mapping modifier names to their constraints, if specified.

"""
function make_sample_dict(samples::NamedTuple, function_specs::NamedTuple)
    constraints = NamedTuple()
    sample_arr = []
    abs_staterror_data = _calc_staterror_data(samples)
    custom = (;)
    for (sample_name, sample) in zip(keys(samples), samples)
        sample_dict =  Dict(key => getfield(sample, key) for key ∈ fieldnames(HistFactorySampleSpec))
        sample_dict[:name] = sample_name   
        mod_arr, temp_constraints, temp_custom = make_modifier_dict(sample, abs_staterror_data, function_specs)
        custom = merge(custom, temp_custom)
        constraints = merge(constraints, temp_constraints)
        sample_dict[:modifiers] = mod_arr
        sample_arr  = push!(sample_arr, sample_dict)
    end
    sample_arr, constraints, custom
end

"""
    _bin_wise_modifier_constraints(constraints::NamedTuple, constraint_parameters)

Merge specific constraints from `constraints` into a new `NamedTuple` based on `constraint_parameters`.

# Arguments
- `constraints`: NamedTuple containing the constraints.
- `constraint_parameters`: Vector of all channel parameters.

# Returns
- New NamedTuple with merged constraints.

# Details
This function matches the constraint parameters against a regular expression pattern to identify specific constraints with the format
"_binN", where N is a number, specifically staterror and shapesys -typed modifiers.
It then merges these constraints into a new NamedTuple using the parameter names as keys.

# Examples
```julia
constraints = (; constraint1 = 10, constraint2 = 20)
constraint_parameters = ["constraint1_bin1", "constraint1_bin2", "constraint2"]

new_constraints = _bin_wise_modifier_constraints(constraints, constraint_parameters)

The resulting new_constraints will be a NamedTuple with the merged constraints: (; constraint1 = 10, constraint2 = 20, constraint1_bin1 = 5, constraint1_bin2 = 15).
"""
function _bin_wise_modifier_contraints(constraints::NamedTuple, constraint_parameters)
    pattern = r"_bin\d+$"  
    nt = (;)
    for element in constraint_parameters
        if match(pattern, string(element)) !== nothing
            name = Symbol(replace(string(element), pattern => s""))
            (constraints[name] !== nothing) ? nt = merge(nt, (Symbol(element) => constraints[name],)) : nothing
        end
    end
    return merge(constraints, nt)
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

#@with_kw struct Custom{T} <: LiteHF.AbstractModifierSpec
#    expression::T
#    function Custom(interp::T) where T 
#        new{T}(interp)
#    end
#end

function make_histfact(histfact_spec::HistFactorySpec, channel_name::Symbol, function_specs=(;))
    sample_arr, constraints, custom = make_sample_dict(histfact_spec.samples, function_specs)
    channel_dict = LiteHF.build_channel(Dict(:samples => sample_arr); misc=Dict())
    modifier_names = unique!(reduce(vcat, [channel_dict[i].modifier_names for i in keys(channel_dict)]))
    channel = LiteHF.build_pyhfchannel((channel_name => channel_dict), modifier_names)
    
    @assert length(modifier_names) == length(keys(channel[2])) #sanity check
    #println("Here", constraints)
    #if constraints != (;)
    #    println("he ", constraints)
    constraints = _bin_wise_modifier_contraints(constraints, keys(channel[2]))
    #end
    constraint_terms = NamedTuple()
    for name in keys(channel[2])
        if (name in keys(constraints)) && (typeof(constraints[name]) != Val{:Const}) 
            constraint_terms = merge(constraint_terms, [name => _make_constraint(constraints[name], channel[2][name])])
        else
            #constraint_terms = merge(constraint_terms, [name => LiteHF._prior(channel[2][name])])
        end
    end
    @info constraint_terms
    HistfactPDF(channel, constraint_terms, modifier_names, custom)
end

