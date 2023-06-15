# This file is a part of HS3.jl, licensed under the MIT License (MIT).
"""
    abstract type AbstractFunctionSpec <: AbstractFunctionalSpec

An abstract type representing a function specification within the HS3 framework.
"""
abstract type AbstractFunctionSpec <: AbstractFunctionalSpec end


"""
    struct FunctionSpec{functype, P<:NamedTuple} <: AbstractFunctionSpec

A type representing a specific function specification within the HS3 framework.

# Fields
- `params::P`: A named tuple containing the parameters specific to the function specification.

# Type Parameters
- `functype`: The specific function type.
- `P<:NamedTuple`: The type of the `params` field.

# Examples
```julia
spec = FunctionSpec{:sum}((summands = [1, 2, 3],))
"""
struct FunctionSpec{functype, P<:NamedTuple} <: AbstractFunctionSpec
    params::P
end

FunctionSpec{functype}(params::P) where {functype,P<:NamedTuple} = FunctionSpec{functype,P}(params)


"""
    parse_and_generate_functionspec(dict::AbstractDict)

Parse a dictionary and delegate the generation of a function specification based on the provided data to the corresponding `generate_functionspec` function.

# Arguments
- `dict::AbstractDict`: A dictionary containing the necessary parameters for the function specification. 
The dictionary should include a `:type` field specifying the function type, along with other relevant fields specific to the function.

# Returns
A function specification of type `FunctionSpec` generated based on the provided data.

# Examples
```julia
dict = Dict(:type => :sum, :summands => [1, 2, 3])
spec = parse_and_generate_functionspec(dict)
"""
function parse_and_generate_functionspec(dict::AbstractDict)
    tp = Val(Symbol(dict[:type]))
    content_dict = filter(entry -> entry[1] != :type, dict)
    generate_functionspec(tp, _typed_content(content_dict))
end


"""
    generate_functionspec(tp::Val, data::NamedTuple)

Generate a function specification of type `FunctionSpec` based on the provided type and data.

# Arguments
- `tp::Val`: A `Val` type representing the specific function type.
- `data::NamedTuple`: A named tuple containing the necessary data for the function specification.

# Returns
A `FunctionSpec` object representing the function specification.

# Examples
```julia
data = (summands = [1, 2, 3, 4, 5])
spec = generate_functionspec(Val(:sum), data)
"""
generate_functionspec(::Val{:sum}, data::NamedTuple) = FunctionSpec{:sum}((summands = data.summands,))

generate_functionspec(::Val{:prod}, data::NamedTuple) = FunctionSpec{:prod}((factors = data.factors,))

# in order for the generic function to work with the topological sorting algorithm the variables need to be extracted
generate_functionspec(::Val{:generic_function}, data::NamedTuple) = FunctionSpec{:generic_function}((expression = Meta.parse(String(_val_content(data.expression))), var = _typed_content([_find_variables_in_expression(Meta.parse(String(_val_content(data.expression))))...]), ))


"""
    generate_functions_specs(spec_array::AbstractArray)

Generate function specifications from an array of function data.

# Arguments
- `spec_array::AbstractArray`: An array containing function data.

# Returns
- `specs::NamedTuple`: A named tuple containing the generated function specifications.

# Examples
```julia
func_data = [
    Dict(:name => "sum_func", :type => "sum", :summands => [1, 2, 3]),
    Dict(:name => "prod_func", :type => "prod", :factors => [4, 5, 6]),
]
func_specs = generate_functions_specs(func_data)

"""
function generate_functions_specs(spec_array::AbstractArray)
    specs = NamedTuple{}()
    for func_data in spec_array
        name = Symbol(func_data.name)
        specs = merge(specs, (name => parse_and_generate_functionspec(func_data),))
    end
    specs
end