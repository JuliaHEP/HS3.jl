# This file is a part of HS3.jl, licensed under the MIT License (MIT).
"""
    make_analyses(ana_spec::AnalysesSpec, specs::NamedTuple)

Create an analyses object based on the provided analysis specification and specifications.

# Arguments
- `ana_spec::AnalysesSpec`: The analysis specification that defines the configuration of the analyses.
- `specs::NamedTuple`: The specifications that contain the required specifiactions of data, functionals etc.

# Returns
- `analyses::NamedTuple`: The analyses object containing likelihood, parameter domain, parameter initialization, prior, and parameters of interest.

"""
function make_analyses(ana_spec::AnalysesSpec, specs::NamedTuple)
    ll_nt = make_likelihood(specs[:likelihoods][Symbol(ana_spec.likelihood)], merge(specs[:functions], specs[:distributions]) ,specs[:data])
    analyses = (
        likelihood = ll_nt.likelihood,
        free_parameters = ll_nt.free_parameters,
        parameter_domain = make_domain_ranges(specs[:domains][Symbol(ana_spec.parameter_domain)],)
    )
    if ana_spec.parameter_init !== nothing
        analyses = merge(analyses, (parameter_init = make_parameterpoints(specs[:parameter_points][Symbol(ana_spec.parameter_init)]),))
    end
    if ana_spec.prior === nothing
        analyses = merge(analyses, (prior = generate_uniform_prior_from_domain(specs[:domains][Symbol(ana_spec.parameter_domain)]),))
    else 
        analyses = merge(analyses, (prior = (specs[:distributions][Symbol(ana_spec.prior)]),))  
    end
    if ana_spec.parameters_of_interest !== nothing
        analyses = merge(analyses, (parameters_of_interest = ana_spec.parameters_of_interest,))
    end
    return analyses;
end
