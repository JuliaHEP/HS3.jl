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

function generate_function(spec::FunctionSpec{:product}) 
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

function generate_function(spec::FunctionSpec{:interpolation0d})
    @assert length(spec.params.codes) == length(spec.params.low) == length(spec.params.high) == length(spec.params.vars)
    f = function (vars)
        total = copy(spec.params.nom)
        for i in 1:length(spec.params.codes)
            total = interpolation0d!(spec.params.codes[i], spec.params.low[i], spec.params.high[i], spec.params.nom, vars[i], total)
        end
        total
    end
    f(spec.params.vars)
end

"""
Implementation of the interpolation0d function, as also used in HistFactory. 
It's depending on an interpolation code, and high, low and nominal values; this is a direct translation of the RooFit code;
a more 'Julia' like implementation would be good.
"""

# As defined in https://root.cern.ch/doc/master/EvaluateFuncs_8h_source.html#l00112 and https://root.cern.ch/doc/master/FlexibleInterpVar_8cxx_source.html#l00196
function interpolation0d!(code, low, high, nominal, paramVal, total)
    if code == 0 
        #piece-wise linear
        if (paramVal > 0)
           return total + paramVal * (high - nominal)
        else
           return total + paramVal * (nominal - low)
        end
    elseif (code == 1) 
        #piece-wise log
        if (paramVal >= 0)
           return total * (high / nominal)^(paramVal)
        else
           return total * (low / nominal)^(-paramVal)
        end
    elseif (code == 2) 
        #parabolic with linear
        a = 0.5 * (high + low) - nominal
        b = 0.5 * (high - low)
        c = 0
        if (paramVal > 1) 
           return total + (2 * a + b) * (paramVal - 1) + high - nominal
        elseif (paramVal < -1) 
           return total + -1 * (2 * a - b) * (paramVal + 1) + low - nominal
        else 
           return total + a * paramVal^2 + b * paramVal + c
        end
    elseif (code == 3) 
      #parabolic version of log-normal
        a = 0.5 * (high + low) - nominal
        b = 0.5 * (high - low)
        c = 0;
        if (paramVal > 1) 
           return total + (2 * a + b) * (paramVal - 1) + high - nominal
        elseif (paramVal < -1) 
           return total + -1 * (2 * a - b) * (paramVal + 1) + low - nominal
        else 
           return total + a * paramVal^2 + b * paramVal + c
        end
    elseif (code == 4) 
        x = paramVal
        boundary = 1.0
        if (x >= boundary) 
           return total * (high / nominal)^(paramVal)
        elseif (x <= -boundary) 
           return total * (low / nominal)^(-paramVal)
        end
        return total * exp(_interpolate_sixt_degree(x, log(low), log(high), log(nominal), boundary))
    else 
        @error "InterpolationCode $code not found/supported"
    end
end

"""
`_interpolate_sixt_degree` 

support function for interpolation of some sixth degree polynome. Also used for HistFactory interpolation.
"""

function _interpolate_sixt_degree(x, low, high, nominal, boundary=1.0)
    t = x / boundary;
    eps_plus =  high - nominal;
    eps_minus = nominal - low;
    S = 0.5 * (eps_plus + eps_minus);
    A = 0.0625 * (eps_plus - eps_minus);
    return x * (S + t * A * (15 + t * t * (-10 + t * t* 3)))
end