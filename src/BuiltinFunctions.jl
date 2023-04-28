# This file is a part of HS3.jl, licensed under the MIT License (MIT).

abstract type AbstractFunctionsSpec <: AbstractFunctionalSpec end

struct FunctionSpec{functype, P<:NamedTuple} <: AbstractFunctionsSpec
    params::P
end

FunctionSpec{functype}(params::P) where {functype,P<:NamedTuple} = FunctionSpec{functype,P}(params)

function AbstractFunctionsSpec(dict::AbstractDict)
    tp = Val(Symbol(dict[:type]))
    content_dict = filter(entry -> entry[1] != :type, dict)
    generate_funcspec(tp, _typed_content(content_dict))
end


#_typed_content(v::Vector{Val}) = Vector(v)
#_typed_content(A::AbstractArray) = [ _typed_content(A[element]) for element in range(1, length(A))]

function _find_variables_in_expression(ex::Expr, vars::Set{Symbol} =Set{Symbol}() )
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


generate_funcspec(::Val{:sum}, data::NamedTuple) = FunctionSpec{:sum}((summands = data.summands,))
function generate_function(spec::FunctionSpec{:sum}) 
    sum(spec.params.summands) 
end

generate_funcspec(::Val{:prod}, data::NamedTuple) = FunctionSpec{:prod}((factors = data.factors,))
function generate_function(spec::FunctionSpec{:prod}) 
    print(spec.params.factors)
    prod(spec.params.factors) 
end

generate_funcspec(::Val{:generic_function}, data::NamedTuple) = FunctionSpec{:generic_function}((expression = Meta.parse(String(_val_content(data.expression))), var = _typed_content([_find_variables_in_expression(Meta.parse(String(_val_content(data.expression))))...]), ))
function generate_function(spec::FunctionSpec{:generic_function}) 
    v = spec.params.var 
    key_names = [_find_variables_in_expression(spec.params.expression)...]
    for i in 1:length(v)
        symbolname = Symbol("$(key_names[i])")
        eval(quote $symbolname = $v[$(i)] end)
    end
    eval(spec.params.expression)
end

@generated function subst_vars(spec::FunctionSpec, v::NamedTuple)
    functype = spec.parameters[1] # e.g. prod
    func_params = spec.parameters[2] # e.g. NamedTuple{(:factors), Tuple{Vector{Real}}}
    pv_type_pairs = map(=>, func_params.parameters[1], func_params.parameters[2].parameters) # e.g. Pair{Symbol, DataType}[:factors => Vector{Real}]
    parnames = func_params.parameters[1] # e.g. (:factors,)
    #remapped_par_content = map(pv -> pv[2] <: Tuple ? map(x -> x <: Val ? :(v.$(x.parameters[1])) : :(spec.params.$(pv[1])), (pv[2].parameters)) : ((pv[2]) <: Val ? :(v.$(pv[2].parameters[1])) : :(spec.params.$(pv[1]))), pv_type_pairs)
    remapped_par_content = Any[]
    parname_iterator = 0
    for pv in pv_type_pairs
        parname_iterator += 1
        #println("hmm", pv[2]<: Vector)
        if pv[2] <: Tuple
            iterator = 0
            help = Vector{Expr}()
            for x in pv[2].parameters 
                iterator +=1                
                if x <: Val
                    push!(help, (:(v.$(x.parameters[1])))) 
                else 
                    push!(help, (:(spec.params.$(pv[1])[$iterator])))
                end
            end 
            push!(remapped_par_content, help)
            #print(iterator)
        else
            if pv[2] <: Val
                push!(remapped_par_content, [:(v.$(pv[2].parameters[1]))]) 
            else
                push!(remapped_par_content, [:(spec.params.$(pv[1]))])
            end
        end
    end
    #println("rem: ", remapped_par_content[1])
    #println("rem: ", remapped_par_content[2])
    #println(typeof(remapped_par_content))

    
    #for pv in pv_type_pairs
    #    #println(pv[2].parameters[1])
    #    if pv[2] <: AbstractArray
    #        pv[2] = Vector(pv[2])
    #        print(pv[2])
    #        for x in pv[2].parameters
    #            #println(x)
    #        end
    #    end
    #end
    #remapped_par_content = 1
    
    #remapped_par_content = map(pv -> pv[2] <: Val ? :(v.$(pv[2].parameters[1].VectorElements )) : :(spec.params.$(pv[1])), pv_type_pairs)
    #remapped_par_content = map(pv -> :(v.(_val_content($(pv[2].parameters[1]))), v.(_val_content($(pv[2].parameters[2])))), pv_type_pairs)length($(pv[2].parameters))
    #remapped_par_content = map(pv -> :((v.$(uffnames))), pv_type_pairs)
    #remapped_par_content = remapped_par_content[1]
    #x = :help
    #print(functype)
    #println(spec.parameters[2].parameters)
    #println("e", remapped_par_content[2])
    #println(pv_type_pairs)
    #remapped_par_content = 1
    tpl_expr = :()
    tpl_expr.args = remapped_par_content[1]
    #tpl_expr2 = :(())    
    #tpl_expr2.args = remapped_par_content[2]
    #print(tpl_expr, tpl_expr2)
    #test = vcat(tpl_expr, tpl_expr2)
    #print(test)
    #println("Args", Meta.show_sexpr(remapped_par_content))
    #tpl_expr.args = (remapped_par_content[1], remapped_par_content[2])
    #tpl_expr.args = :(())
    #tpl_expr.args.args = remapped_par_content
    #for i = 1:length(remapped_par_content)
    
    #println("w", length(remapped_par_content[2]))
    #for element in remapped_par_content
    #    print(collect(element...))
    #    push!(tpl_expr.args, collect(element...))
    #end
    #try 
    #    tpl_expr.args = remapped_par_content[1] 
    #catch 
    #    try 
    #        tpl_expr.args = remapped_par_content
    #        print("aaaaa")
    #    catch
    #        throw("The Function is not readable")
    #    end
    #end 
    ##print(dump(tpl_expr))
    #new_params = :(NamedTuple{$(QuoteNode(parnames[1]))}(($tpl_expr)))
    #new_params2 = :(NamedTuple{$(QuoteNode(parnames[2]))}(($tpl_expr2)))
    #print(tpl_expr.args[1])
    #new_params = :[($(QuoteNode(parnames[1])) => $tpl_expr)]
    if length(remapped_par_content[1]) == 1
        new_params = :[($(QuoteNode(parnames[1])) => $tpl_expr[1])]
    else 
        new_params = :[($(QuoteNode(parnames[1])) => $tpl_expr)]
    end
    #new_params2 = :[($(QuoteNode(parnames[2])) => $tpl_expr2[1])]
    #print(new_params)
    #print(new_params2)
    #nt = NamedTuple()
    h = :((; $new_params...))
    #h = :(merge($h, $new_params2))
    #K = :($h)
    #print(K)

    #print(b)
    #b
    l = length(remapped_par_content)
    for var_num = 2:l
        tpl_expr = :()
        tpl_expr.args = remapped_par_content[var_num]
        if length(remapped_par_content[var_num]) == 1
            new_params = :[($(QuoteNode(parnames[var_num])) => $tpl_expr[1])]
        else 
            new_params = :[($(QuoteNode(parnames[var_num])) => $tpl_expr)]
        end
        h = :(merge($h, $new_params))
        #print("dljlsaökjdj", tpl_expr)
    end
    b = :(FunctionSpec{$(QuoteNode(functype))}($h))
    b
end

function make_funcspecs(func_array::AbstractArray)
    specs = NamedTuple()
    for func in 1:length(func_array)
        specs = merge(specs, (Symbol(func_array[func].name) => AbstractFunctionsSpec(func_array[func]),))
    end
    :($specs)
end


make_functional(spec::AbstractFunctionsSpec) = generate_function ∘ Base.Fix1(subst_vars, spec)
#make_func(spec::FunctionSpec{:generic_function}) = Base.Fix1(generate_function, spec)

export make_func
export AbstractFunctionsSpec
export subst_func_vars
export make_funcspecs