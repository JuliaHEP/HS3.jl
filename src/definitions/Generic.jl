struct Generic{T<:Real} <: ContinuousUnivariateDistribution
    formula::Expr
    parameters::Vector{Any}
    var::Symbol
    parameter_names::Vector{Symbol}
    func::Function
    #a::T  # lower bound of the interval
    #b::T  # upper bound of the interval
end
#
## Constructor for Generic from a string
function remove_val_and_corresponding_symbol(array_with_val, array_with_symbols)
    val_index = findfirst(x -> typeof(x).name.wrapper == Val, array_with_val)
    if val_index !== nothing
        deleteat!(array_with_val, val_index)
        deleteat!(array_with_symbols, val_index)
    end
    return array_with_val, array_with_symbols
end

#function Generic(formula_expr::Expr, parameters::Vector{Any}, var::Symbol)  
#    parameter_names = [_find_variables_in_expression(formula_expr)...]
#    remove_val_and_corresponding_symbol(parameters, parameter_names)
#    print(var)
#    func = x->begin
#        for i in 1:length(parameters)
#            symbolname = Symbol("$(parameter_names[i])")
#            eval(quote $symbolname = $parameters[$(i)] end)
#        end
#        #println(eval(formula_expr))
#        symbolname = Symbol("$(var)")
#        eval(quote $symbolname = x end)
#        eval(formula_expr)
#    end
#    return Generic{Float64}(formula_expr, parameters, var, parameter_names, func)
#end
#
function _substitute(expr::Expr, from::Symbol, to::Symbol)
    for i in eachindex(expr.args)
        if isa(expr.args[i], Symbol)
            if expr.args[i] == from
                expr.args[i] = to
            end
        elseif isa(expr.args[i], Expr)
            _substitute(expr.args[i], from, to)
        end
    end
    return expr
end

function Generic(formula_expr::Expr, parameters::Vector{Any}, var::Symbol)
    parameter_names = collect(_find_variables_in_expression(formula_expr))
    remove_val_and_corresponding_symbol(parameters, parameter_names)
    assignments = Expr(:block)
    for (name, value) in zip(parameter_names, parameters)
        push!(assignments.args, :( $name = $value ))
    end

    func_expr = Expr(:block, assignments, formula_expr)
    func_expr = _substitute(func_expr, var, :x)
    func = eval(:(x -> $func_expr))

    return Generic{Float64}(formula_expr, parameters, var, parameter_names, func)
end

function pdf(d::Generic, x::Real) 
    @eval($d.func($x))
end

logpdf(d::Generic, x::Real) = log(pdf(d, x))
#
## Convert the stored expression into a callable function
#function make_function(d::Generic, var::Symbol)
#    @variables x
#    subs_expr = _substitute(d.formula, Dict(var => x))
#    parameters = [x; [k for k in keys(d.parameters)]]
#    func_expr = :($(subs_expr))
#    func = eval(:(($parameters...)-> $func_expr))
#    return (v) -> func(v, values(d.parameters)...)
#end
#
## Implement the PDF using the constructed function
#function pdf(d::Generic{T}, x::T) where T
#    key_names = [_find_variables_in_expression(spec.params.expression)...]
#    for i in 1:length(d.parameters)
#        symbolname = Symbol("$(key_names[i])")
#        eval(quote $symbolname = $v[$(i)] end)
#    end
#    eval(spec.params.expression)
#    #if x < d.a || x > d.b
#    #    return 0.0
#    #else
#    func = y -> 
#    println(func)
#    return func(x)
#    #end
#end
##
## Implement the log PDF using the constructed function
#function logpdf(d::Generic{T}, var::Symbol, x::T) where T
#    #if x < d.a || x > d.b
#    #    return -Inf
#    #else
#    func = make_function(d, var)
#    return log(func(x))
#    #end
#end
#struct Generic <: ContinuousUnivariateDistribution
#    expression::Expr
#    parameters::Set{Symbol}
#    var::Symbol
#end
#
##function Generic(expression::Expr, parameters::Set{Symbol}; check_args::Bool=true)
##    @check_args Generic (parameters, length(parameters) >= 1)
##    return Generic(expression, parameters)
##end
#
#Generic(expression::Expr, domain::AbstractArray{<:AbstractRange}, variates::Set{Symbol}; check_args::Bool=true) = Generic(expression,_find_variables_in_expression(expression), domain, variates, _n_val(Val{:Generic}(), variates, Meta.parse(expression), domain))
#Generic(expression::String, domain::AbstractArray{<:AbstractRange}, variates::Set{Symbol}; check_args::Bool=true) = Generic(Meta.parse(expression), _find_variables_in_expression(Meta.parse(expression)), domain, variates, _n_val(Val{:Generic}(), variates, Meta.parse(expression), domain))
##Generic(expression::Expr, min::Number, max::Number; check_args::Bool=true) = Generic(expression,_find_variables_in_expression(expression), (min:max))
##Generic(expression::String, min::Number, max::Number; check_args::Bool=true) = Generic(Meta.parse(expression), _find_variables_in_expression(Meta.parse(expression)), (min:max))
##Polynomial(a::AbstractArray{Real}, min::Number, max::Number; check_args::Bool=true) = Polynomial(a, (min:max); check_args=check_args)
#
###Polynomial(a::AbstractArray{Real}; check_args::Bool=true) = Polynomial(a, (0:1); check_args=check_args)
##Polynomial() = Polynomial([1., 1., 1.], (0:1); check_args=false)
#
#params(d::Generic) = (d.expression,)
#
#expression(d::Generic) = d.expression
#
#
#parameters(d::Generic) = d.parameters
#
#
#function _not_norm_pdf(d::Generic, x::Vector)
#    for (i, name) in enumerate(d.variates)
#        #(Meta.parse("$name = :(x[$i])" ))
#        (eval(Meta.parse("$name = $(x[i])")))
#    end
#    #println(d.expression)
#    result = eval(d.expression)
#    return result
#end
#
#function _pdf_to_vegas(expression, variates, x::Vector)
#    for (i, name) in enumerate(variates)
#        #(Meta.parse("$name = :(x[$i])" ))
#        (eval(Meta.parse("$name = $(x[i])")))
#    end
#    #println(d.expression)
#    result = eval(expression)
#    return result
#end
#
#pdf(d::Generic, x::Vector) = _not_norm_pdf(d, x) /d._norm
#logpdf(d::Generic, x::Vector) = log(pdf(d, x))
#
#function _make_integrand(expression, variates, bounds, dims::Integer)
#    return function integrand(x, f)
#        f[1] = _transform_to_fixed_boundaries( x -> _pdf_to_vegas(expression, variates, x), dims, first.(bounds), last.(bounds))(x)
#    end 
#end 
#
#function _transform_to_fixed_boundaries(f, dim, a::Vector, b::Vector)
#    return x -> begin
#        jacobian = 1.0
#        transformed_x = Vector{Float64}(undef, dim)
#        
#        for i in 1:dim
#            transformed_x[i] = a[i] + (b[i] - a[i]) * x[i]
#            jacobian *= (b[i] - a[i])
#        end
#        return f(transformed_x) * jacobian
#    end
#end
#
#function _n_val(::Val{:Generic}, variates, expression, bounds)
#    dims = length(variates)
#    integrand = _make_integrand(expression, variates, bounds, dims)
#    vegas(integrand, dims, 1, atol=1e-3, rtol=1e-3)[1][1]
#end
##function normalize_expression(expr_str, variable_limits)
##    # Define symbols dynamically based on the number of variables
##    var_symbols = [Symbol(var[1]) for var in variable_limits]
##    @variables vars[var_symbols...]
##    
##    # Convert the expression string into a symbolic expression
##    expr = Symbolics.parse_expression(expr_str)
##    
##    # _substitute the variables in the expression
##    expr_sub = Symbolics._substitute(expr, vars)
##    
##    # Define the limits
##    limits = [(vars[symbol], limit[2], limit[3]) for (symbol, limit) in zip(var_symbols, variable_limits)]
##
##    # Compute the integral over the entire space
##    total_mass = foldl((acc, limit) -> integral(acc, limit...), expr_sub, init = limits)
##
##    # Normalize the expression
##    normalized_expr = expr_sub / total_mass
##
##    return normalized_expr
##end
##
#