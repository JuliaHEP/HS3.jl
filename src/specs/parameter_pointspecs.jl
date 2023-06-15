# This file is a part of HS3.jl, licensed under the MIT License (MIT).
"""
    struct ParameterPointSpec <: AbstractParameterPointsSpec
A specification for a single parameter point in a parameter space.

# Fields
- `value`: The numerical value of the parameter.
- `isconst`: A boolean flag indicating if the parameter is constant (default: `false`).

"""
@with_kw struct ParameterPointSpec <: AbstractHS3Spec
    value::Number
    isconst::Bool = false 
end


"""
    generate_parameterpoint_spec(nt)

Generate a ParameterPointSpec object from the given named tuple.

# Arguments
- `nt::NamedTuple`: The named tuple containing the parameter point information.
  - `value`: The value of the parameter.
  - `const` (optional): The const flag indicating if the parameter is constant.

# Returns
A `ParameterPointSpec` object representing the parameter point.

"""
generate_parameterpoint_spec(nt::NamedTuple{(:value,), <:Tuple{Number}}) = ParameterPointSpec(value = nt.value, isconst = false) 

generate_parameterpoint_spec(nt::NamedTuple{(:value, :const), <:Tuple{Number, Integer}}) = ParameterPointSpec(value = nt.value, isconst = nt[:const]) 

"""
    parse_set_of_parameter_points(arr::AbstractArray)
Generate parameter point specifications based on the provided array.

# Arguments
- `arr`: Abstract array containing the parameter point data.

# Returns
A named tuple containing the generated parameter point specifications.

"""
function parse_set_of_parameter_points(arr::AbstractArray)
    specs = NamedTuple{}()
    for element in arr
        temp = NamedTuple(filter(entry -> entry[1] != :name, element))
        specs = merge(specs, (Symbol(element.name) => generate_parameterpoint_spec(temp),))
    end
    return specs
end


"""
    generate_parameter_point_specs(param_array::AbstractArray)
Generate parameter point specifications by generating parameter points from the provided array of parameter point sets.

# Arguments
- `param_array`: Abstract array containing the parameter sets.

# Returns
A named tuple containing the generated parameter point specifications.

"""
function generate_parameter_points_specs(param_array::AbstractArray)
    specs = NamedTuple()
    for param_set in param_array
        specs = merge(specs, (Symbol(param_set.name) => parse_set_of_parameter_points(param_set.parameters, ),))
    end
    specs
end