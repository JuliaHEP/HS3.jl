# This file is a part of HS3.jl, licensed under the MIT License (MIT).

### For now this is just a dump of code
function HS3_output()
    open("my_new_file.json", "w") do io
        JSON3.pretty(io, a)
    end
end

for expr in code
    if expr isa Expr && expr.head == :call
        operation = push!(operation, expr.args[1].name)
        for arg in expr.args
            if typeof(arg) == QuoteNode
                params_names = push!(param_names, arg.value)
            end
        end
    end
end

code = (code_lowered(logdensityof(a))[1].code)