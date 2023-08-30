# This file is a part of HS3.jl, licensed under the MIT License (MIT).
"""
    substitute_variables(spec::AbstractFunctionalSpec, v::NamedTuple)

Substitute variables in the given functional specification with provided values.

# Arguments
- `spec::AbstractFunctionalSpec`: The functional specification.
- `v::NamedTuple`: The variables/parameters to be substitutet with values.

# Returns
The modified functional specification with substituted variables.

"""
@generated function substitute_variables(spec::AbstractFunctionalSpec, v::NamedTuple)
    functional_type = spec.parameters[1] # e.g. poisson
    #varnames = spec.parameters[2] # e.g. (:x,), name of the variate
    functional_parameters = spec.parameters[2]# e.g. NamedTuple{(:lambda,), Tuple{Val{:param_lambda}}}
    parameter_value_pairs = map(=>, functional_parameters.parameters[1], functional_parameters.parameters[2].parameters) # e.g. Pair{Symbol, DataType}[:lambda => Val{:param_lambda}]
    parameter_names = functional_parameters.parameters[1] # e.g. (:lambda,)
    num_parameters = length(parameter_names)
    remapped_parameter_content = Any[]
    for pv in parameter_value_pairs
        if pv[2] <: Tuple #check if parameter is array-like
            iterator = 0
            vector_expr = Vector{Expr}()
            for x in pv[2].parameters #iterate over all parameters in array
                iterator +=1                
                if x <: Val && x.parameters[1] in [fieldnames(v)...] #check if parameter is already defined or needs to be replaced with element of v
                    push!(vector_expr, (:(v.$(x.parameters[1])))) 
                else 
                    push!(vector_expr, (:(spec.params.$(pv[1])[$iterator])))
                end
            end 
            push!(remapped_parameter_content, vector_expr)
        else
            if pv[2] <: Val && pv[2].parameters[1] in [fieldnames(v)...] #check if parameter is already defined or needs to be replaced with element of v
                push!(remapped_parameter_content, [:(v.$(pv[2].parameters[1]))]) 
            else
                push!(remapped_parameter_content, [:(spec.params.$(pv[1]))])
            end
        end
    end
    tpl_expr = :()
    tpl_expr.args = remapped_parameter_content[1] # e.g. :(v.param_lambda)
    if length(remapped_parameter_content[1]) == 1
        new_params = :[($(QuoteNode(parameter_names[1])) => $tpl_expr[1])]
    else 
        new_params = :[($(QuoteNode(parameter_names[1])) => $tpl_expr)]
    end
    merged_params = :((; $new_params...))
    num_parameters = length(remapped_parameter_content)
    for var_num = 2:num_parameters
        tpl_expr = :()
        tpl_expr.args = remapped_parameter_content[var_num]
        if length(remapped_parameter_content[var_num]) == 1
            new_params = :[($(QuoteNode(parameter_names[var_num])) => $tpl_expr[1])]
        else 
            new_params = :[($(QuoteNode(parameter_names[var_num])) => $tpl_expr)]
        end
        merged_params = :(merge($merged_params, $new_params))
    end
    if spec <: DistributionSpec
        :(DistributionSpec{$(QuoteNode(functional_type))}($merged_params))
    elseif spec <: FunctionSpec
       :(FunctionSpec{$(QuoteNode(functional_type))}($merged_params))
    end
end


"""
    make_functional(spec, sorted_functionals=[], level=0)

Recursively create a functional representation of the given specification, optionally under consideration of dependencies of other functionals.

# Arguments
- `spec`: The functional specification.
- `sorted_functionals`: (Optional) The sorted functionals used in the creation of the functional representation.
- `level`: (Optional) The level of the specification in the functional hierarchy.

# Returns
An executeable functional representation of the specification.

"""

function make_functional(spec::AbstractFunctionalSpec, sorted_functionals::Vector, level::Int64 = 0)
    return params -> begin
        functionals = NamedTuple()
        sorted_functionals = filter(x -> x[1] > level, sorted_functionals)
        for funct in reverse(sorted_functionals) 
            functionals = merge(functionals, (first(last(funct)) => make_functional(funct[2][2], sorted_functionals, level+1)(params), ))
        end
        #println("hereeee: ", functionals)
        funct = make_functional(spec)
        #println(zip((first.(last.(sorted_functionals))), functionals))
        params = merge(params, zip((first.(last.(sorted_functionals))), functionals))
        #println("bla ", sorted_functionals)
        #for (k, v) in zip(keys(functionals), functionals)
        #    println(functionals)
        #    (typeof(v) <: Number) && break
        #    functionals = merge(functionals, (k => v(params),))
        #end
        #@info (funct, " what ", level, " ", funct(merge(functionals, params)))
        funct((params))
    end
end

function make_functional(spec::Union{FunctionSpec{:generic_function}} , sorted_functionals::Vector, level::Int64=0)
    functionals = NamedTuple()
    vars = spec.params.var
    sorted_functionals = filter(x -> x[1] > level, sorted_functionals)
    for funct in reverse(sorted_functionals) 
        if funct in vars 
            functionals = merge(functionals, (first(last(funct)) => make_functional(funct[2][2]), ))
        end
    end
    funct = make_functional(spec)
    return params -> begin
        for (k, v) in zip(keys(functionals), functionals)
            (typeof(v) <: Number) && break
            functionals = merge(functionals, (k => v(merge(functionals, params)),))
        end
        funct(merge(functionals, params))
    end
end


make_functional(spec::AbstractDistributionSpec) = generate_distribution ∘ Base.Fix1(substitute_variables, spec)

make_functional(spec::AbstractFunctionSpec) = generate_function ∘ Base.Fix1(substitute_variables, spec)

#Redundant function, use only for validation purposes!
make_functional(spec::HistFactorySpec) = make_histfact(spec, :default_name)