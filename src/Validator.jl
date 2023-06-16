# This file is a part of HS3.jl, licensed under the MIT License (MIT).
"""
validate_file(path::String)

Validate a file containing specifications.

# Arguments
- `path`: A string specifying the path to the file.

# Returns
- Nothing

"""
function validate_file(path::String)
    #check if JSON-File is valid
    dict = file_to_dict(path)
    #check correct implementation of specs
    specs = generate_specs(dict)
    functionals = (;)
    if haskey(specs, :analyses)
        for analysis in specs[:analyses]
            make_analyses(analysis, specs)
        end
        for likelihood in specs[:likelihoods]
            if haskey(specs, :functions)
                functionals = merge(specs[:functions], specs[:distributions])
                make_likelihood(likelihood, functionals, specs[:data])
            else
                functionals = specs[:distributions]
                make_likelihood(likelihood, functionals, specs[:data])
            end
        end
        for data in specs[:data]
            make_data(data)
        end
        sorted_functionals = topological_sort(functionals)
        for functional in functionals
            make_functional(functional, sorted_functionals)
        end
    end
end

export validate_file
