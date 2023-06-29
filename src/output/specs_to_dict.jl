# This file is a part of HS3.jl, licensed under the MIT License (MIT).

function get_variable_name(var::Any)
    for (name, value) in Base.iterators(vars(Main))
        if value === var
            return string(name)
        end
    end
    error("HS3 variable not found")
end

_type_to_HS3(::Type{<:HS3.BinnedDataSpec{N}}) where N = "binned"
_type_to_HS3(::Type{HS3.HistFactorySpec}) = "histfactory_dist"
_type_to_HS3(::Type{HS3.NormsysSpec}) = "normsys"
_type_to_HS3(::Type{HS3.StaterrorSpec}) = "staterror"
_type_to_HS3(::Type{HS3.NormfactorSpec}) = "normfactor"



function specs_to_json_dict(myStruct::NormsysSpec, name::String = "default")
    jsonDict = Dict()   
    jsonDict["data"] = Dict("hi" => getfield(myStruct, :data).hi, "lo" => getfield(myStruct, :data).lo)
    jsonDict["type"] = "normsys"
    getfield(myStruct, :constraint) !== nothing ? jsonDict["constraint"] = getfield(myStruct, :constraint) : nothing 
    jsonDict["name"] = name
    jsonDict
end

function specs_to_json_dict(myStruct::HistFactorySampleSpec, name::String="default")
    jsonDict = Dict()   
    if getfield(myStruct, :errors) !== nothing 
        jsonDict["data"] =  Dict("contents" => getfield(myStruct, :data), "errors" => getfield(myStruct, :errors))
    else
        jsonDict["data"] =  Dict("contents" => getfield(myStruct, :data))
    end
    jsonDict["modifiers"] = specs_to_json_dict(getfield(myStruct, :modifiers))
    jsonDict["name"] = name
    return jsonDict
end

function specs_to_json_dict(myStruct::AbstractHS3Spec, name::String="default")
    fields = fieldnames(typeof(myStruct))
    jsonDict = Dict()   
    for field in fields
        if getfield(myStruct, Symbol(field)) !== nothing
            if typeof(getfield(myStruct, Symbol(field))) <: Union{Number, String}
                jsonDict[String(field)] = getfield(myStruct, Symbol((field)))
            else
                jsonDict[String(field)] = specs_to_json_dict(getfield(myStruct, Symbol(field)))
            end
        end
    end
    if any(x -> isa(myStruct, x), [AbstractFunctionalSpec, AbstractModifierSpec, AbstractDataSpec])
        jsonDict["type"] = _type_to_HS3(typeof(myStruct))
    end
    
    jsonDict["name"] = _val_content(name)
    return jsonDict
end

function specs_to_json_dict(val::Any, name::String)
    if typeof(name) <: Val 
        name = _val_content(name)
    end
    if typeof(val) <: Val
        val = _val_content(val)
    end
    return name => val
end

function specs_to_json_dict(val::Any)
    if typeof(val) <: Val
        val = _val_content(val)
    end
    return val
end

function specs_to_json_dict(nt::NamedTuple)
    return [specs_to_json_dict((v), String(k)) for (k, v) in zip(keys(nt), nt)]
end

function specs_to_json_dict(nt::NamedTuple, toplevel::Bool)
    spec_types = [
        :analyses,
        :likelihoods,
        :distributions,
        :parameter_points,
        :functions,
        :data,
        :domains,
        :metadata,
        :misc]
    jsonDict = Dict()
    if toplevel == false 
        specs_to_json_dict(nt::NamedTuple)
    else
        if all(x -> x ∈ spec_types, keys(nt))
            for (k, v) in zip(keys(nt), nt)
                jsonDict[String(k)] = specs_to_json_dict(v)
            end
            return jsonDict
        else
            @error "The provided Named Tuple is not a top-level HS3 representation"
        end
    end
end

#function specs_to_json_dict(nt::NamedTuple)
#    jsonDict = Dict()
#    spec_types = [
#        :analyses,
#        :likelihoods,
#        :distributions,
#        :parameter_points,
#        :functions,
#        :data,
#        :domains,
#        :metadata,
#        :misc]
#    end
#    if keys(nt) .∈ specs_types
#        specs_to_json_dict()
#    else
#end