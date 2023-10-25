# This file is a part of HS3.jl, licensed under the MIT License (MIT).
"""
    struct HistfactPDF

The `HistfactPDF` struct represents a probability density function (PDF) defined within the context of a HistFactory model given by LiteHF.jl.

# Fields
- `channel::Tuple`: A tuple representing the channel information associated with the PDF.
- `prior::NamedTuple`: A named tuple representing the prior/modifier information associated with the PDF.

"""
struct HistfactPDF 
    channel::Function 
    prior::NamedTuple 
    order::Vector
    customs::NamedTuple
end
#
#
#"""
#    make_dict_from_modifier(modifier) -> Dict
#
#Create a dictionary representation of a given `modifier` object usable for the LiteHF.jl Package.
#
## Arguments
#- `modifier`: The input object.
#
## Returns
#A dictionary representing the input `modifier` object.
#
#"""
#make_dict_from_modifier(modifier::NormfactorSpec, sample::HistFactorySampleSpec, abs_staterror_data::Tuple{Vector{Float64}, Vector{Float64}}) = Dict(:type => "normfactor", :data => nothing)
#
#make_dict_from_modifier(modifier::NormsysSpec, sample::HistFactorySampleSpec, abs_staterror_data::Tuple{Vector{Float64}, Vector{Float64}}) = Dict(:type => "normsys", :data => Dict(:hi => modifier.data.hi, :lo => modifier.data.lo))
#
#function make_dict_from_modifier(modifier::StaterrorSpec, sample::HistFactorySampleSpec, abs_staterror_data::Tuple{Vector{Float64}, Vector{Float64}}) 
#    @assert length(abs_staterror_data[2]) == length(abs_staterror_data[1])
#    ##@info abs_staterror_data
#    ##abs_staterror_data = 100.
#    #@info sample.errors
#    #@info abs_staterror_data
#    σ = copy(abs_staterror_data[1])
#    for i in 1:length(abs_staterror_data[1]) 
#        if abs_staterror_data[2][i] != 0
#            σ[i] = abs_staterror_data[1][i] / abs_staterror_data[2][i]
#        else 
#            σ[i] = sqrt(abs_staterror_data[1][i])
#        end 
#    end
#    Dict(:type => "staterror", :data =>  σ)
#end
#make_dict_from_modifier(modifier::HistosysSpec, sample::HistFactorySampleSpec, abs_staterror_data::Tuple{Vector{Float64}, Vector{Float64}}) = Dict(:type => "histosys", :data => Dict(:hi_data => modifier.data.hi.contents, :lo_data => modifier.data.lo.contents))
#
#function make_dict_from_modifier(modifier::ShapesysSpec, sample::HistFactorySampleSpec, abs_staterror_data::Tuple{Vector{Float64}, Vector{Float64}})
#    if modifier.data === nothing 
#        data = zeros(length(sample.data))
#    else 
#        data = modifier.data.vals 
#    end
#    Dict(:type => "shapesys", :data => data)
#end 
#
#make_dict_from_modifier(modifier::ShapefactorSpec, sample::HistFactorySampleSpec, abs_staterror_data::Tuple{Vector{Float64}, Vector{Float64}}) = Dict(:type => "shapefactor", :data => modifier.data)
#
#make_dict_from_modifier(modifier::CustomModifierSpec, sample::HistFactorySampleSpec, abs_staterror_data::Tuple{Vector{Float64}, Vector{Float64}}) = Dict(:type => "custom", :data => nothing)

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

#
#function make_modifier_dict(sample::HistFactorySampleSpec, abs_staterror_data::Tuple{Vector{Float64}, Vector{Float64}}, function_specs::NamedTuple)
#    constraints = NamedTuple()
#    mod_arr = []
#    custom = (;)
#    sorted = topological_sort(function_specs)
#    custom_parameter_names = []
#    custom_mod_name = [] #TODO here better stuff ends to go
#    for (mod_name, mod) in zip(keys(sample.modifiers), sample.modifiers)
#        for element in mod 
#            #@info element
#            if typeof(element) != CustomModifierSpec 
#                mod_dict = make_dict_from_modifier(element, sample, abs_staterror_data)
#                mod_dict[:name] = String(mod_name)
#                mod_arr = push!(mod_arr, mod_dict)
#            else
#                func = make_functional(function_specs[Symbol(mod_name)], sorted)
#                @info func
#                #TODO: Account for non generic functions
#                parameter_names = []
#                if typeof(function_specs[Symbol(mod_name)]) == FunctionSpec{:generic_function}
#                    parameter_names = _val_content(collect(function_specs[Symbol(mod_name)].params.var))
#                elseif typeof(function_specs[Symbol(mod_name)]) <: FunctionSpec{:interpolation0d}
#                    #println("yes ", parameter_names)
#                    parameter_names = _val_content(collect(function_specs[Symbol(mod_name)].params.vars))
#                    #@info (parameter_names)
#                    #sleep(10)
#                end
#                #@info parameter_names
#
#                #println(parameter_names)
#                for i in function_specs 
#                    if typeof(i) == FunctionSpec{:generic_function}
#                        parameter_names = append!(parameter_names, _val_content(collect(i.params.var)))
#                    elseif typeof(i) == FunctionSpec{:interpolation0d}
#                        parameter_names = append!(parameter_names, _val_content(collect(i.params.vars)))
#                    end
#                end
#                #@info parameter_names
#                mod_dict = Dict(:type=> "custom", :data =>(params -> func(params), parameter_names) , :name => mod_name)
#                mod_arr = push!(mod_arr, mod_dict)
#                custom_parameter_names = push!(custom_parameter_names, parameter_names)
#                custom_mod_name = push!(custom_mod_name, mod_name)
#                #custom = merge(custom, (mod_name => parameter_names,))
#            end
#            if :constraint ∈ fieldnames(typeof(element)) && element.constraint !== nothing
#                #constraint_names = push!(constraint_names, mod_dict[name])
#                constraints = merge(constraints, [Symbol(mod_dict[:name]) => element.constraint,])
#            end 
#            if :constraint_name ∈ fieldnames(typeof(element)) && element.constraint_name !== nothing
#                #constraint_names = push!(constraint_names, mod_dict[name])
#                #@info "----------------------------------------HEEERREEEEE -------------------------------"
#                constraints = merge(constraints, [Symbol(mod_dict[:name]) => function_specs[Symbol(_val_content(element.constraint_name))],])
#            end 
#        end
#    end
#    #@info constraints
#    mod_arr, constraints, (; zip(custom_mod_name, custom_parameter_names)...)
#end


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

#_make_constraint(::Val{:Gauss}, modifier::Normsys) =  Distributions.Normal(0, 1.)
#
#_make_constraint(::Val{:Gauss}, modifier::Histosys) = Distributions.Normal(0, 1.)
#
#function _make_constraint(::Val{:Gauss}, modifier::Shapesys) 
#    #@info "Shapesys: " sqrt(modifier.σn2) modifier.moddata
#    #@info (modifier.σn2)
#    #@info (modifier.moddata)
#    Distributions.Normal(0 , 1/sqrt(modifier.moddata))
#end
#
#_make_constraint(::Val{:Gauss}, modifier::Staterror) = Distributions.Normal(0, modifier.σ)
#
#_make_constraint(::Val{:Poisson}, modifier::Normsys) = RelaxedPoisson(1.)
#
#_make_constraint(::Val{:Poisson}, modifier::Histosys) = RelaxedPoisson(1.)
#
#_make_constraint(::Val{:Poisson}, modifier::Shapesys) = RelaxedPoisson(Float64(modifier.moddata))
#
#_make_constraint(::Val{:Poisson}, modifier::Staterror) = RelaxedPoisson(Float64(modifier.σ^2))

#_make_constraint(::Val{:Const}, modifier::AbstractModifier) = FlatPrior(1, 1) # TODO: this needs improvement

#_make_constraint(::Val{:Const}, modifier::Shapesys) = RelaxedPoisson(0.)# TODO: this needs improvement
#_make_constraint(dist::Any, modfier::Normfactor) = Distributions.truncated(Distributions.Normal(1, dist.params.sigma), 0, 10)

#_make_constraint(::Any, ::Any) = @error "Not a valid entry and thus not defined"

#_make_constraint(::Val{:LogNormal} ) = Distributions.LogNormal()


#function _calc_staterror_data(samples::NamedTuple)
#    data = zero(samples[1].data)
#    errors = zero(samples[1].data)
#    for sample in samples
#        for mods in values(sample.modifiers)
#            if any(mod -> typeof(mod) == StaterrorSpec, mods)
#                data .+= (sample.data)
#                errors .+= (sample.errors).^2 
#            end
#        end
#    end
#    @assert length(data) == length(errors)
#    data, sqrt.(errors)
#end

#function merge_constraints(nt1::NamedTuple, nt2::NamedTuple, mod_arr::AbstractArray)
#    pairs_to_merge = Pair{Symbol, Any}[]  # An array to hold pairs we want to merge
#
#    # Start with the entries in nt1
#    for (k, v) in pairs(nt1)
#        push!(pairs_to_merge, k => v)
#    end
#
#    # Then consider the entries in nt2
#    for (k, v) in pairs(nt2)
#        if haskey(nt1, k) 
#            # If the key already exists in nt1 and its value is Val{:P}
#            if nt1[k] == Val{:Poisson}() && v == Val{:Const}()
#                # Skip this pair
#                mod_arr = filter(dict -> dict[:name] != String(k), mod_arr)
#                continue
#            end
#        end
#        # Either add a new key or update an existing one
#        idx = findfirst(p -> p.first == k, pairs_to_merge)
#        if idx === nothing
#            push!(pairs_to_merge, k => v)
#        else
#            pairs_to_merge[idx] = k => v
#        end
#    end
#    return NamedTupleTools.namedtuple(pairs_to_merge), mod_arr
# end

"""
make_sample_dict(samples)

Create a sample dictionary and constraint dictionary from a collection of samples.

# Arguments
- `samples`: A collection of samples.

# Returns
- `sample_arr`: An array of sample dictionaries.
- `constraints`: A dictionary mapping modifier names to their constraints, if specified.

"""
#function make_sample_dict(samples::NamedTuple, function_specs::NamedTuple)
#    constraints = NamedTuple()
#    sample_arr = []
#    abs_staterror_data = _calc_staterror_data(samples)
#    custom = (;)
#    for (sample_name, sample) in zip(keys(samples), samples)
#        sample_dict =  Dict(key => getfield(sample, key) for key ∈ fieldnames(HistFactorySampleSpec))
#        sample_dict[:name] = sample_name   
#        mod_arr, temp_constraints, temp_custom = make_modifier_dict(sample, abs_staterror_data, function_specs)
#        #println(temp_constraints)
#        custom = merge(custom, temp_custom)
#        constraints, mod_arr = merge_constraints(constraints, temp_constraints, mod_arr)
#        sample_dict[:modifiers] = mod_arr
#        sample_arr  = push!(sample_arr, sample_dict)
#    end
#    sample_arr, constraints, custom
#end

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
#function _bin_wise_modifier_contraints(constraints::NamedTuple, constraint_parameters)
#    pattern = r"_bin_\d+$"  
#    nt = (;)
#    for element in constraint_parameters
#        if match(pattern, string(element)) !== nothing
#            println(element)
#            name = Symbol(replace(string(element), pattern => s""))
#            #println(constraints)
#            (constraints[name] !== nothing) ? nt = merge(nt, (Symbol(element) => constraints[name],)) : nothing
#        end
#    end
#    #println(constraints, constraint_parameters)
#    return merge(constraints, nt)
#end


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

function _create_custom_nt(customs::AbstractArray)
    if length(customs) > 0
        list = [(first(item)) for item in customs]
    #@info NamedTupleTools.namedtuple(list)
        return NamedTupleTools.namedtuple(list)
    else
        return (;)
    end
    #values_list = [(values(item)) for item in list]
    #println(keys_list, values_list)
    #return NamedTuple{tuple(keys_list...)}(values_list)
end


function make_histfact(histfact_spec::HistFactorySpec, name::Symbol, function_specs=(;))
    #sample_arr, constraints, custom = make_sample_dict(histfact_spec.samples, function_specs)
    #println(constraints)
    channel_dict  = build_channel(histfact_spec, name, function_specs)
    #modifier_names = unique!(reduce(vcat, [channel_dict[i].modifier_names for i in keys(channel_dict)]))
    channel = _create_expected(channel_dict)
    #sleep(20)
    #println(dump(channel[1]))
    #@assert length(modifier_names) == length(keys(channel[2])) #sanity check
    #println("Here", constraints)
    #if constraints != (;)
    #    println("he ", constraints)
    #@info keys(channel[2])
    #constraints = _bin_wise_modifier_contraints(constraints, keys(channel[2]))
    #println("there", constraints)
    #end
    #constraint_terms = NamedTuple()
    #for name in keys(channel[2])
    #    if (name in keys(constraints)) && (typeof(constraints[name]) != Val{:Const}) 
    #        constraint_terms = merge(constraint_terms, [name => _make_constraint(constraints[name], channel[2][name])])
    #    end
    #end
    #println(constraint_terms)
    #@info (channel_dict[4])
    #sleep(30)
    custom = _create_custom_nt(channel_dict[4])
    order = setdiff(channel[2], collect(keys(custom)))
    #@info keys(custom)
    HistfactPDF(channel[1], channel[3], order, custom)
end

