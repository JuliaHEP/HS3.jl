abstract type AbstractParameterPointsSpec  <: AbstractHS3Spec end

function generate_parameterpoint_specs(arr::AbstractArray)
    nt = NamedTuple()
    for element in arr
        temp = NamedTuple(filter(entry -> entry[1] != :name, element))
        #temp = NamedTuple(filter(entry -> entry[1] != :name, temp))
        nt = merge(nt, (Symbol(element.name) => gen_point_specs(temp),),)
    end
    return nt
end

Parameters.@with_kw struct ParameterPointSpec <: AbstractParameterPointsSpec 
    value::Number
    isconst::Bool = false 
end

#ParameterPointSpec(value::Number, constant::Integer) = ParameterPointSpec(value, Bool(constant))
gen_point_specs(nt::NamedTuple{(:value,), <:Tuple{Number}}) = ParameterPointSpec(value = nt.value, isconst = false) 
gen_point_specs(nt::NamedTuple{(:value, :const), <:Tuple{Number, Integer}}) = ParameterPointSpec(value = nt.value, isconst = nt[:const]) 
#
function make_parameterpointspecs(param_array::AbstractArray)
    specs = NamedTuple()
    for param_set in param_array
        specs = merge(specs, (Symbol(param_set.name) => generate_parameterpoint_specs(param_set.parameters, ),))
    end
    specs
end

function make_parameterpoints(pointspecs::NamedTuple)
    nt = NamedTuple()
    for (k,v) in zip(keys(pointspecs), pointspecs)
        nt = merge(nt, (k => v.value,))
    end
    nt
end