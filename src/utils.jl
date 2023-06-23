# This file is a part of HS3.jl, licensed under the MIT License (MIT).

"""
    _find_variables_in_expression(ex::Expr, vars::Set{Symbol} = Set{Symbol}())

Recursively searches for variables in an expression.

# Arguments
- `ex::Expr`: The expression to search for variables.
- `vars::Set{Symbol}`: Set of variables found so far (default: empty set).

# Returns
- `vars::Set{Symbol}`: Set of variables found in the expression, excluding common operators (+, -, *, /).

"""
function _find_variables_in_expression(ex::Expr, vars::Set{Symbol} = Set{Symbol}() )
    for arg in ex.args
        if isa(arg, Symbol)
            # If the argument is a symbol, add it to the set of variables
            push!(vars, arg)
        elseif isa(arg, Expr)
            # If the argument is another expression, recursively search for variables
            _find_variables_in_expression(arg, vars)
        end
    end
    return filter(x -> x ∉ [:+, :-, :*, :/], vars)
end


"""
    _typed_content(data)

Converts the input data to a typed representation.

# Arguments
- `data`: The input data to convert.

# Returns
- Typed representation of the input data.

The `_typed_content` function applies specific transformations based on the type of the input `data`:
- If `data` is a `String` or a `Symbol`, it is converted to a `Val` wrapping the corresponding `Symbol`.
- If `data` is a `Number`, it is returned as is.
- If `data` is a collection of objects, e.g. an array or a dict, it recursively applies `_typed_content` to each element and returns the resulting collection.
- For any other input, it returns `nothing`.
"""

function _typed_content(dict::AbstractDict)
    data = collect(dict)
    NamedTuple{(map(first, data)...,)}(map(_typed_content ∘ last, data))
end

_typed_content(str::String) = Val(Symbol(str))

_typed_content(s::Symbol) = Val(s)

_typed_content(x::Number) = x

_typed_content(A::AbstractArray{<:Number}) = A

_typed_content(A::AbstractArray) = (_typed_content.(A)...,)

_typed_content(::Nothing) = nothing

_typed_content(A::NamedTuple) = A

"""
    _val_content

Extracts the content from a `Val` or an `AbstractArray`.

# Arguments
- `::Val{x}`: A `Val` object containing a value `x` *or*
- `A::AbstractArray`: An `AbstractArray` object.

# Returns
- If the input is a `Val`, it returns the value `x`.
- If the input is an `AbstractArray`, it recursively applies `_val_content` to each element and returns the resulting array.
- For any other input, it returns `nothing`.

"""

_val_content(::Val{x}) where x = x

_val_content(A::AbstractArray) = _val_content.(A)

_val_content(A::Any) = A