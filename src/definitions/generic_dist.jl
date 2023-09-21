"""
Generic Distribution 
Provides Distribution with generic pdf definition through expression
"""

using Distributions: @check_args, @distr_support, ContinuousMultivariateDistribution
struct Generic <: ContinuousMultivariateDistribution
    expression::Expr
    parameters::Set{Symbol}
    domain::AbstractArray{AbstractRange}
    variates::Set{Symbol}
    _norm::Float64
end

#function Generic(expression::Expr, parameters::Set{Symbol}; check_args::Bool=true)
#    @check_args Generic (parameters, length(parameters) >= 1)
#    return Generic(expression, parameters)
#end

Generic(expression::Expr, domain::AbstractArray{<:AbstractRange}, variates::Set{Symbol}; check_args::Bool=true) = Generic(expression,_find_variables_in_expression(expression), domain, variates, _n_val(Val{:Generic}(), variates, Meta.parse(expression), domain))
Generic(expression::String, domain::AbstractArray{<:AbstractRange}, variates::Set{Symbol}; check_args::Bool=true) = Generic(Meta.parse(expression), _find_variables_in_expression(Meta.parse(expression)), domain, variates, _n_val(Val{:Generic}(), variates, Meta.parse(expression), domain))
#Generic(expression::Expr, min::Number, max::Number; check_args::Bool=true) = Generic(expression,_find_variables_in_expression(expression), (min:max))
#Generic(expression::String, min::Number, max::Number; check_args::Bool=true) = Generic(Meta.parse(expression), _find_variables_in_expression(Meta.parse(expression)), (min:max))
#Polynomial(a::AbstractArray{Real}, min::Number, max::Number; check_args::Bool=true) = Polynomial(a, (min:max); check_args=check_args)

##Polynomial(a::AbstractArray{Real}; check_args::Bool=true) = Polynomial(a, (0:1); check_args=check_args)
#Polynomial() = Polynomial([1., 1., 1.], (0:1); check_args=false)

params(d::Generic) = (d.expression,)

expression(d::Generic) = d.expression


parameters(d::Generic) = d.parameters


function _not_norm_pdf(d::Generic, x::Vector)
    for (i, name) in enumerate(d.variates)
        #(Meta.parse("$name = :(x[$i])" ))
        (eval(Meta.parse("$name = $(x[i])")))
    end
    #println(d.expression)
    result = eval(d.expression)
    return result
end

function _pdf_to_vegas(expression, variates, x::Vector)
    for (i, name) in enumerate(variates)
        #(Meta.parse("$name = :(x[$i])" ))
        (eval(Meta.parse("$name = $(x[i])")))
    end
    #println(d.expression)
    result = eval(expression)
    return result
end

pdf(d::Generic, x::Vector) = _not_norm_pdf(d, x) /d._norm
logpdf(d::Generic, x::Vector) = log(pdf(d, x))

function _make_integrand(expression, variates, bounds, dims::Integer)
    return function integrand(x, f)
        f[1] = _transform_to_fixed_boundaries( x -> _pdf_to_vegas(expression, variates, x), dims, first.(bounds), last.(bounds))(x)
    end 
end 

function _transform_to_fixed_boundaries(f, dim, a::Vector, b::Vector)
    return x -> begin
        jacobian = 1.0
        transformed_x = Vector{Float64}(undef, dim)
        
        for i in 1:dim
            transformed_x[i] = a[i] + (b[i] - a[i]) * x[i]
            jacobian *= (b[i] - a[i])
        end
        return f(transformed_x) * jacobian
    end
end

function _n_val(::Val{:Generic}, variates, expression, bounds)
    dims = length(variates)
    integrand = _make_integrand(expression, variates, bounds, dims)
    vegas(integrand, dims, 1, atol=1e-3, rtol=1e-3)[1][1]
end
#function normalize_expression(expr_str, variable_limits)
#    # Define symbols dynamically based on the number of variables
#    var_symbols = [Symbol(var[1]) for var in variable_limits]
#    @variables vars[var_symbols...]
#    
#    # Convert the expression string into a symbolic expression
#    expr = Symbolics.parse_expression(expr_str)
#    
#    # Substitute the variables in the expression
#    expr_sub = Symbolics.substitute(expr, vars)
#    
#    # Define the limits
#    limits = [(vars[symbol], limit[2], limit[3]) for (symbol, limit) in zip(var_symbols, variable_limits)]
#
#    # Compute the integral over the entire space
#    total_mass = foldl((acc, limit) -> integral(acc, limit...), expr_sub, init = limits)
#
#    # Normalize the expression
#    normalized_expr = expr_sub / total_mass
#
#    return normalized_expr
#end
#
