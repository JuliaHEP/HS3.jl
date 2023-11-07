"""
`_calc_channel_data`:
Support function needed for 'Staterror' modifier; calculates the complete data and associated error in a channel
    """  

function _calc_channel_data(samples::NamedTuple)
    data = zero(samples[1].data)
    errors = zero(samples[1].data)
    for sample in samples
        for mods in values(sample.modifiers)
            if any(mod -> typeof(mod) == StaterrorSpec, mods)
                corrected_errors = ifelse.(sample.errors .== 0, sqrt.(abs.(sample.data)), sample.errors)
                data .+= sample.data
                errors .+= corrected_errors.^2
            end
        end
    end
    @assert length(data) == length(errors)
    data, sqrt.(errors)
end


""" 
`build_channel`
Main function to create HistFactory typed objects
    """

function build_channel(histfact_spec, channel_name, function_specs)
    expectedcounts = []
    names = []
    constraints = []
    channel_data = _calc_channel_data(histfact_spec.samples)
    customs = []
    for (sample_name, sample) in zip(keys(histfact_spec.samples), histfact_spec.samples)
        e, n, c, custom = build_sample(sample, function_specs, channel_data, channel_name) #dict?, is sample_name relevant?
        push!(expectedcounts, e)
        push!(names, n)
        push!(constraints, c)
        push!(customs, custom)
    end
    customs = filter(x -> !(x isa Array && isempty(x)), customs)
    expectedcounts, names, constraints, customs
end


"""
Support function to assess the staterror modifier names; they are depended on bin contents and other.
    This needs to be included in HS3 documentation!
    """
function _assert_modifier_names!(names, modifier, name, channel_name)
    push!(names, name)
end

function _assert_modifier_names!(names, modifier::Vector{<:Staterror}, name, channel_name)
    channel_name = String(channel_name)
    channel_name = replace(channel_name, "model_" => "")
    for i in eachindex(modifier)
        push!(names, Symbol("gamma_stat_", channel_name, "_bin_$(i-1)")) ##improve this
    end
    return names
end


"""
Building function for a single sample; returns a function that represents the expected yield. 
    Handles the association of modifiers 
    """ 

function build_sample(sample, function_specs, channel_data, channel_name)
    modifier_arr = []
    names = Symbol[]
    constraint_arr = []
    customs = []
    for (name, mod) in zip(keys(sample.modifiers), sample.modifiers)
        modifiers = []
        constraints = []
        for m in mod
            modifier_i, constraint_i, custom_i = build_modifier(m, name; sample.data, sample.errors, function_specs, channel_data)
            push!(modifiers, modifier_i)
            push!(constraints, constraint_i)
            push!(customs, custom_i)
            names = _assert_modifier_names!(names, modifier_i, name, channel_name)
        end
        push!(modifier_arr, modifiers)
        push!(constraint_arr, constraints)
        #_push_if_not_nothing!(constraints, )
    end
    modifier_arr = any(x->x <: Vector, typeof.(modifier_arr)) ? reduce(vcat,  modifier_arr) : modifier_arr #flatten it
    modifier_arr = any(x->x <: Vector, typeof.(modifier_arr)) ? reduce(vcat,  modifier_arr) : modifier_arr
    constraint_arr = any(x->x <: Vector, typeof.(constraint_arr)) ?  reduce(vcat, reduce(vcat, constraint_arr)) : constraint_arr
    customs = filter(x -> x !== nothing, customs)
    ExpCounts(collect(sample.data), names, modifier_arr), names, constraint_arr, customs
end

"""
Create constraints for modifiers.

    """
_make_constraint(::Val{:Gauss}, ::Val{:Histosys}) = Distributions.Normal(0, 1.)

_make_constraint(::Val{:Poisson}, ::Val{:Histosys}) = RelaxedPoisson(1.)

function build_modifier(modobj::HistosysSpec, modname; data, errors, function_specs, channel_data)
    Histosys(data, modobj.data.lo.contents, modobj.data.hi.contents), _make_constraint(modobj.constraint, Val(:Histosys)), nothing
end

_make_constraint(::Val{:Gauss}, ::Val{:Normsys}) =  Distributions.Normal(0, 1.)

_make_constraint(::Val{:Poisson}, ::Val{:Normsys}) = RelaxedPoisson(1.)

function build_modifier(modobj::NormsysSpec, modname; data, errors, function_specs, channel_data)
    Normsys(data, modobj.data.hi, modobj.data.lo), _make_constraint(modobj.constraint, Val(:Normsys)), nothing
end

_make_constraint(::Val{:Gauss}, σ,  ::Val{:Staterror}) = Distributions.Normal(0, σ)

_make_constraint(::Val{:Poisson}, σ, ::Val{:Staterror}) = RelaxedPoisson(Float64(σ^2))

function build_modifier(modobj::StaterrorSpec, modname; data, errors, function_specs, channel_data)
    σ = copy(channel_data[1])
    for i in 1:length(channel_data[1]) 
        if channel_data[2][i] != 0
            σ[i] = channel_data[1][i] / channel_data[2][i]
        else 
            σ[i] = sqrt(channel_data[1][i])
        end 
    end
    Staterror.(σ, length(data), eachindex(σ)),  _make_constraint.(modobj.constraint, σ, Val(:Staterror)), nothing
end

function _make_constraint(::Val{:Gauss}, val, ::Val{:Shapesys}) #maybe also unconstrained?
    Distributions.Normal(0 , 1/sqrt(val))
end

_make_constraint(::Nothing, val, ::Val{:Shapesys}) = nothing

_make_constraint(::Val{:Poisson}, val, ::Val{:Shapesys}) = RelaxedPoisson(Float64(val))

function build_modifier(modobj::ShapesysSpec, modname; data, errors, function_specs, channel_data)
    vals = []
    if modobj.data === nothing 
        vals = zeros(length(sample.data))
    else 
        vals = modifier.data.vals 
    end
    Shapesys.((data ./ vals).^2, vals, nbins, eachindex(vals)), _make_constraint.(modobj.constraint, vals, Val(:Shapesys)), nothing
end

function build_modifier(modobj::CustomModifierSpec, modname; data, errors, function_specs, channel_data)
    CustomMod(params -> make_functional(function_specs[Symbol(modname)])(params), _parameter_names_from_function(function_specs[Symbol(modname)])), nothing, (modname => _parameter_names_from_function(function_specs[Symbol(modname)])) 
end

_make_constraint(constraint_name, function_specs,  ::Val{:Normfactor}) = Distributions.Normal(1, (function_specs[Symbol(_val_content(constraint_name))].params.sigma))

_make_constraint(::Nothing, function_specs, ::Val{:Normfactor}) = nothing

function build_modifier(modobj::NormfactorSpec, modname; data, errors, function_specs, channel_data)
    Normfactor(), _make_constraint(modobj.constraint_name, function_specs, Val(:Normfactor)), nothing
end

_parameter_names_from_function(function_spec::FunctionSpec{:generic_function}) = _val_content(collect(function_spec.params.var))

_parameter_names_from_function(function_spec::FunctionSpec{:interpolation0d}) = _val_content(collect(function_spec.params.vars))

function _create_constraint_nt(names::Vector{Symbol}, entries::Vector)
    @assert length(names) == length(entries)
    valid_pairs = ((n, e) for (n, e) in zip(names, entries) if e !== nothing)
    return NamedTuple(Dict(valid_pairs))
end

function _extract_pairs(arr, pairs=[])
    for item in arr
        if item isa Pair
            push!(pairs, item)
        elseif item isa Array
            _extract_pairs(item, pairs)
        end
    end
    return pairs
end

function _create_expected(channel)
    constraint_arr = []
    (expcounts, all_names, constraints, customs) = channel
    constraints = _create_constraint_nt(reduce(vcat, all_names), reduce(vcat, constraints))
    custom_pairs = _extract_pairs(customs)
    custom_dict = Dict(custom_pairs)
    all_modified_names = [ [get(custom_dict, elem, elem) for elem in names] for names in all_names]
    all_unique_names = unique!(reduce(vcat, reduce(vcat, all_modified_names)))
    masks = Tuple([(typeof(i) == Symbol ? findfirst(==(i), all_unique_names) : [findfirst(==(j), all_unique_names) for j in i]) for i in names] for names in all_modified_names)
    expected = (αs) -> internal_expected(Tuple(expcounts), masks, αs)
    return expected, all_unique_names, constraints
end

struct ExpCounts{M}
    nominal::Vector{Float64}
    modifier_names::Vector{Symbol}
    modifiers::M
end

ExpCounts(nominal::Vector{Float64}, names::Vector{Symbol}, modifiers::AbstractVector) = ExpCounts(nominal, names, tuple(modifiers...))

(E::ExpCounts)(αs::AbstractArray) =_expkernel(E.modifiers, E.nominal, αs)

"""
Calculate the actual expectation
"""
@generated function _expkernel(modifiers::Tuple{Vararg{Any, N}}, nominal::Vector{Float64}, αs::AbstractArray) where N
    code = quote
        additive = zeros(Float64, length(nominal))
        factor = ones(Float64, length(nominal))
    end

    for i = 1:N
        push!(code.args, quote
            @inbounds exp_mod!(modifiers[$i].interp, additive, factor, αs[$i])
        end)
    end

    push!(code.args, :(return factor .* (additive .+ nominal)))

    return code
end

function Base.show(io::IO, ::Type{<:ExpCounts}) 
    println(io, "ExpCounts{}")
end

for T in (Normsys, Normfactor, Histosys, Staterror, Shapefactor, CustomMod)
    @eval function Base.show(io::IO, E::$T)
        interp = Base.typename(typeof(E.interp)).name
        print(io, $T)
        #print(io, $T, "{$interp}")
    end
end

#function Base.show(io::IO, E::ExpCounts)
#    modifiers = E.modifiers
#    #elip = length(modifiers) > 5 ? "...\n" : ""
#    println(io, "\n ExpCounts: Nominal: $(E.nominal),")
#    print.("($(E.modifier_names[i]) = $(E.modifiers[i])), " for i in 1:length(E.modifier_names))
#    println("")
#end


 function exp_mod!(interp::PiecewiseInterpolation, additive::Vector{Float64}, factor::Vector{Float64}, α::Float64)
    additive .+= interp(α)
    return
end

function exp_mod!(interp::Union{InterpCode4, typeof(identity),Function}, additive::Vector{Float64}, factor::Vector{Float64}, α::Float64)
    factor .*= interp(α)
    return
end

#function exp_mod!(interp, additive::Vector{Float64}, factor::Vector{Float64}, α)
#    factor .*= interp(α)
#    return
#end

function exp_mod!(interp::CustomInterpolation, additive::Vector{Float64}, factor::Vector{Float64}, α::Vector{Float64})
    factor .*= interp(α)
    return
end

@generated function internal_expected(Es::Tuple, Vs::Tuple, αs)::Float64
    #@assert Es <: Tuple
    @views expand(i) = i == 1 ? :(Es[1](αs[Vs[1]])) : :(+(Es[$i](αs[Vs[$i]]), $(expand(i-1))))
    return expand(length(Es.parameters))
end

function pyhf_loglikelihoodof(expected::Function, obs::Vector{Float64}) 
    f(e, o) = e < zero(e) ? oftype(e, -Inf) : _relaxedpoislogpdf(e, o)
    return function log_likelihood(αs::AbstractArray)::Float64
        mapreduce(+, expected(αs), obs) do E, O # sum over channels
            sum(Base.Broadcast.broadcasted(f, E, O))
        end
    end
end

@generated function internal_constrainteval(pris::Tuple, αs::Vector{Float64})::Float64
    #@assert pris <: Tuple
    @views expand(i) = i == 1 ? :(logpdf(pris[1], αs[1])) : :(+(logpdf(pris[$i], αs[$i]), $(expand(i-1))))
    return expand(length(pris.parameters))
end

function pyhf_logpriorof(priors)
    pris = values(priors)
    function (αs)
        #@assert length(αs) == length(pris)
        internal_constrainteval(pris, αs)
    end
end