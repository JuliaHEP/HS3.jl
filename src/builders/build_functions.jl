# This file is a part of HS3.jl, licensed under the MIT License (MIT).
"""
    generate_function(spec::FunctionSpec)

Generate a function based on the given function specification.

# Arguments
- `spec::FunctionSpec`: The function specification.

# Returns
The generated function based on the function specification.

"""

function generate_function(spec::FunctionSpec{:sum}) 
    sum(spec.params.summands) 
end

function generate_function(spec::FunctionSpec{:prod}) 
    prod(spec.params.factors) 
end

function generate_function(spec::FunctionSpec{:generic_function}) 
    v = spec.params.var 
    key_names = [_find_variables_in_expression(spec.params.expression)...] # TODO: use params.var
    for i in 1:length(v)
        symbolname = Symbol("$(key_names[i])")
        eval(quote $symbolname = $v[$(i)] end)
    end
    eval(spec.params.expression)
end