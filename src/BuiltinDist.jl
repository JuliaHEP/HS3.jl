# This file is a part of HS3.jl, licensed under the MIT License (MIT).
abstract type AbstractFunctionalSpec <: AbstractHS3Spec end

abstract type AbstractDistributionsSpec <: AbstractFunctionalSpec end

struct DistributionSpec{disttype, P<:NamedTuple, T<:NamedTuple} <: AbstractDistributionsSpec
    params::P
    variables::T
end

DistributionSpec{disttype}(params::P, variables::T) where {disttype, P<:NamedTuple, T<:NamedTuple} = DistributionSpec{disttype,P, T}(params, variables)
DistributionSpec{disttype}(params::P) where {disttype, P<:NamedTuple} = DistributionSpec{disttype,P,typeof((;))}(params,(;))

function _typed_content(dict::Dict)
    data = collect(dict)
    NamedTuple{(map(first, data)...,)}(map(_typed_content ∘ last, data))
end

_typed_content(str::String) = Val(Symbol(str))

_typed_content(s::Symbol) = Val(s)

_typed_content(x::Number) = x

_typed_content(A::AbstractArray{<:Number}) = A

_typed_content(A::AbstractArray) = (_typed_content.(A)...,)

_typed_content(A::JSON3.Object) = _typed_content(Dict(A))

_typed_content(A::PropDicts.PropDict) = _typed_content(Dict(A))


function generate_distspec(dict::AbstractDict) #TODO: Falscher Name
    tp = Val(Symbol(dict[:type]))
    content_dict = filter(entry -> entry[1] != :type, dict)
    generate_distspec(tp, _typed_content(content_dict))
end

_val_content(::Val{x}) where x = x
_val_content(A::AbstractArray) = _val_content.(A)
_val_content(::Any)  = nothing

generate_distspec(::Val{:exponential_dist}, data::NamedTuple) = DistributionSpec{:exponential_dist}((c = data.c,), (x = data.x,))
function generate_distribution(spec::DistributionSpec{:exponential_dist})
    Distributions.Exponential(-1/spec.params.c)
end

generate_distspec(::Val{:normal_dist}, data::NamedTuple) = generate_distspec(Val(:gaussian_dist), data)
generate_distspec(::Val{:gaussian_dist}, data::NamedTuple) = DistributionSpec{:gaussian_dist}((mean = data.mean, sigma = data.sigma,), (x = data.x,))
function generate_distribution(spec::DistributionSpec{:gaussian_dist})
    Distributions.Normal(spec.params.mean, spec.params.sigma)
end

generate_distspec(::Val{:lognormal_dist}, data::NamedTuple) = DistributionSpec{:lognormal_dist}((mean = data.mean, k = data.k,), (x = data.x,))
function generate_distribution(spec::DistributionSpec{:lognormal_dist})
    Distributions.LogNormal(spec.params.mean, spec.params.k)
end

generate_distspec(::Val{:multinormal_dist}, data::NamedTuple) = DistributionSpec{:multinormal_dist}((covariances = reduce(hcat, data.covariances), mean = data.mean,), (x = data.x,) )
function generate_distribution(spec::DistributionSpec{:multinormal_dist}) 
    Distributions.MvNormal(spec.params.mean, spec.params.covariances)
end#

generate_distspec(::Val{:poisson_dist}, data::NamedTuple) = DistributionSpec{:poisson_dist}((mean = data.mean,), (x = data.x,))
function generate_distribution(spec::DistributionSpec{:poisson_dist}) 
    Distributions.Poisson(spec.params.mean)
end

generate_distspec(::Val{:argus_dist}, data::NamedTuple) = DistributionSpec{:argus_dist}((resonance = data.resonance, slope = data.slope, power = data.power,), (x = data.mass,))
function generate_distribution(spec::DistributionSpec{:argus_dist})
    @error "not yet implemented"
end

generate_distspec(::Val{:mixture_dist}, data::NamedTuple) = DistributionSpec{:mixture_dist}((summands = data.summands, coefficients = data.coefficients,)) 
function generate_distribution(spec::DistributionSpec{:mixture_dist}) 
    return Distributions.MixtureModel(collect(spec.params.summands), collect(spec.params.coefficients)/(sum(collect((spec.params.coefficients))))) 
   
        #res = 0.0
        #for (dist, coeff) in zip(spec.params.summands, spec.params.coefficients) 
        #    res += coeff * dist
        #end
        #res

end

generate_distspec(::Val{:product_dist}, data::NamedTuple) = DistributionSpec{:product_dist}((factors = data.factors, )) 
function generate_distribution(spec::DistributionSpec{:product_dist}) 
    return Distributions.product_distribution(collect(spec.params.factors)...)
end

generate_distspec(::Val{:uniform_dist}, data::NamedTuple) = DistributionSpec{:uniform_dist}((min = data.min, max = data.max,))
function generate_distribution(spec::DistributionSpec{:uniform_dist})
    return Distributions.Uniform(spec.params.min, spec.params.max)
end

generate_distspec(::Val{:polynomial_dist}, data::NamedTuple) = DistributionSpec{:polynomial_dist}((coefficients = data.coefficients,), (x= data.x, ))
function generate_distribution(spec::DistributionSpec{:polynomial_dist}) 
    @warn "not yet implemented"
end#




generate_distspec(::Val{:CrystalBall_dist}, data::NamedTuple) = DistributionSpec{:CrystalBall_dist}((alpha = data.alpha, x= data.x, m0 = data.m0, n = data.n, sigma = data.sigma, ))
function generate_distribution(spec::DistributionSpec{:CrystalBall_dist}) 
    @warn "not yet implemented and not yet part of HS3"
end



@generated function subst_vars(spec::DistributionSpec, v)
    disttype = spec.parameters[1] # e.g. poisson
    #varnames = spec.parameters[2] # e.g. (:x,), name of the variate
    dist_params = spec.parameters[2] # e.g. NamedTuple{(:lambda,), Tuple{Val{:param_lambda}}}
    #### do we need the variate here? spec.parameters[3])
    pv_type_pairs = map(=>, dist_params.parameters[1], dist_params.parameters[2].parameters) # e.g. Pair{Symbol, DataType}[:lambda => Val{:param_lambda}]
    ## dist_params.parameters[1] e.g. (:lambda,)
    ## dist_params.parameters[2].parameters e.g. svec(Val{:param_lambda})
    parnames = dist_params.parameters[1] # e.g. (:lambda,)
    #remapped_par_content = map(pv -> pv[2] <: Val ? :(v.$(pv[2].parameters[1])) : :(spec.params.$(pv[1])), pv_type_pairs) ## e.g. Expr[:(v.param_lambda)]
    ##spec.params.$(pv[1]) : e.g. (:lambda => 2)
    remapped_par_content = Any[]
    parname_iterator = 0
    for pv in pv_type_pairs
        parname_iterator += 1
        if pv[2] <: Tuple #check if parameter is array-like
            iterator = 0
            vect = Vector{Expr}()
            for x in pv[2].parameters 
                iterator +=1                
                if x <: Val #check if parameter is already defined or needs to be replaced with element of v
                    push!(vect, (:(v.$(x.parameters[1])))) 
                else 
                    push!(vect, (:(spec.params.$(pv[1])[$iterator])))
                end
            end 
            push!(remapped_par_content, vect)
        else
            if pv[2] <: Val #check if parameter is already defined or needs to be replaced with element of v
                push!(remapped_par_content, [:(v.$(pv[2].parameters[1]))]) 
            else
                push!(remapped_par_content, [:(spec.params.$(pv[1]))])
            end
        end
    end
    tpl_expr = :()
    tpl_expr.args = remapped_par_content[1] # e.g. :(v.param_lambda)
    if length(remapped_par_content[1]) == 1
        new_params = :[($(QuoteNode(parnames[1])) => $tpl_expr[1])]
    else 
        new_params = :[($(QuoteNode(parnames[1])) => $tpl_expr)]
    end
    h = :((; $new_params...))
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
    end
    b = :(DistributionSpec{$(QuoteNode(disttype))}($h))
    println(b)
    b
end

function make_functional(spec::AbstractFunctionalSpec, sorted_dists::Vector, level::Int64 = 0)
    return params -> begin
        functionals = []
        sorted_dists = filter((x) -> x[1] > level, sorted_dists)
        for dist in sorted_dists 
            functionals = push!(functionals, make_functional(dist[2][2], sorted_dists, dist[1])(params))
        end
        params = merge(params, zip((first.(last.(sorted_dists))), functionals))
        new_spec = subst_vars(spec, params)
        make_functional(new_spec)(params)
    end 
end

make_functional(spec::AbstractDistributionsSpec) = generate_distribution ∘ Base.Fix1(subst_vars, spec)
#
#function make_markov_kernel(spec::DistributionSpec{:mixture_dist}, dists::NamedTuple)
#    return v -> begin
#        markov_kernels = []
#        for dist in dists 
#            markov_kernels = push!(markov_kernels, make_markov_kernel(dist)(v))
#        end
#        v = merge(v, zip(keys(dists), markov_kernels))
#        spec = subst_vars(spec, v)
#        make_markov_kernel(spec)(v)
#    end
#end
#
#function make_markov_kernel(spec::DistributionSpec{:product_dist}, dists::NamedTuple)
#    return v -> begin
#        markov_kernels = []
#        for dist in dists 
#            markov_kernels = push!(markov_kernels, make_markov_kernel(dist)(v))
#        end
#        v = merge(v, zip(keys(dists), markov_kernels))
#        spec = subst_vars(spec, v)
#        make_markov_kernel(spec)(v)
#    end
#end

function make_distspecs(dist_array::JSON3.Array)
    specs = NamedTuple()
    for dist in 1:length(dist_array)
        specs = merge(specs, (Symbol(dist_array[dist].name) => generate_distspec(dist_array[dist]),))
    end
    :($specs)
end

export make_distspecs

#function generate_prior_from_distribution(spec::AbstractDistributionsSpec, marginal_specs::NamedTuple)
#    for dist in spec.params.factors
#        _val_content(marginal_specs[Symbol(dist)].variables.x)  
#    end
#end