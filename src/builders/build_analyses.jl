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
    if haskey(specs, :functions)
        functional_specs = merge(specs[:functions], specs[:distributions])
    else
        functional_specs = specs[:distributions]
    end
    ll_nt = make_likelihood(specs[:likelihoods][Symbol(ana_spec.likelihood)], functional_specs, specs[:data])
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

"""
make_analyses(path_to_file::String, name_of_analysis::String)

Create an analysis from a file and return the result.

# Arguments
- `path_to_file`: A string specifying the path to the file.
- `name_of_analysis`: A string specifying the name of the analysis.

# Returns
- The analysis result based on the specified file and analysis name.

"""
function make_analyses(path_to_file::String, name_of_analysis::String)
    specs = generate_specs(path_to_file) #maybe global storage?
    return make_analyses(specs.analyses[Symbol(name_of_analysis)], specs)
end
