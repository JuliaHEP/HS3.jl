# This file is a part of HS3.jl, licensed under the MIT License (MIT).
abstract type AbstractAnalysesSpec <: AbstractHS3Spec end

Parameters.@with_kw struct AnalysesSpec <:AbstractAnalysesSpec
    likelihood::String 
    parameter_domain::String 
    parameters_of_interest::Union{AbstractArray{String}, Nothing} = nothing
    parameter_init::Union{String, Nothing} = nothing
    prior::Union{String, Nothing} = nothing
end


function make_analysesspecs(analyses_array::AbstractArray)
    specs = NamedTuple()
    for ana in analyses_array
        specs = merge(specs, (Symbol(ana.name) => AnalysesSpec(; NamedTuple(filter(entry -> entry[1] != :name, ana))...),),)
    end
    specs
end

function make_analyses(ana_spec::AnalysesSpec, specs::NamedTuple)
    analyses = (likelihood = make_likelihood(specs[:likelihoods][Symbol(ana_spec.likelihood)], specs[:functionals], specs[:data]),
                parameter_domain = make_domain(specs[:domains][Symbol(ana_spec.parameter_domain)], ))
    if ana_spec.parameter_init !== nothing
        analyses = merge(analyses, (parameter_init = make_parameterpoints(specs[:parameter_points][Symbol(ana_spec.parameter_init)]),))
    end
    if ana_spec.prior === nothing
        analyses = merge(analyses, (prior = generate_prior_from_domain(specs[:domains][Symbol(ana_spec.parameter_domain)]),))
    else 
        @warn "not yet implemented"
        #TODO: check how prior would be implemented in HS3.
        #analyses = merge(analyses, (prior = generate_prior_from_distribution(specs[:funcitonals][Symbol(ana_spec.parameter_domain)]),))  
    end
    if ana_spec.parameters_of_interest !== nothing
        analyses = merge(analyses, (parameters_of_interest = ana_spec.parameters_of_interest,))
    end
    return analyses
end

_make_specs(::Val{:metadata}, dict) = return dict
_make_specs(::Val{:domains}, dict) = generate_domainspecs(dict)


_make_specs(::Any, dict) = nothing 
function generate_all_specs(d::AbstractDict)
    nt = NamedTuple()
    for k in keys(d)
        nt = merge(nt, (k => _make_specs(Val(k), d[k]),))
    end
    nt
end

function generate_all_specs(path_to_file::String)
    generate_all_specs(file_to_dict(path_to_file))
end