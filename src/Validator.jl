# This file is a part of HS3.jl, licensed under the MIT License (MIT).

function validate_file(path::String)
    #check if JSON-File is valid
    dict = file_to_dict(path)
    #check correct implementation of specs
    specs = generate_specs(dict)
    if haskey(specs, :analyses)
        for analysis in specs[:analyses]
            make_analyses(analysis, specs)
        end
        for likelihood in specs[:likelihoods]
            make_likelihood(likelihood, merge(specs[:functions], specs[:distributions]), specs[:data])
        end
        for data in specs[:data]
            make_data(data)
        end
        sorted_functionals = topological_sort(merge(specs[:distributions], specs[:functions]))
        for functional in merge(specs[:distributions], specs[:functions])
            make_functional(functional, sorted_functionals)
        end
    end
end

export validate_file
