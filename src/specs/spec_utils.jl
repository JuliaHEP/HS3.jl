# This file is a part of HS3.jl, licensed under the MIT License (MIT).
"""
generate_specs(spec_dict)

Generate the specifications for different components of HS3 based on the provided dictionary.

# Arguments
- `spec_dict::AbstractDict`: The dictionary containing the specifications for different HS3 components.
Alternatively the path to the HS3 file might be directly passed to this function. 

# Returns
A named tuple containing the generated specifications for each component.

"""
function generate_specs(spec_dict::AbstractDict)
    spec_types = [
        :analyses,
        :likelihoods,
        :distributions,
        :parameter_points,
        :functions,
        :data,
        :domains,
        :metadata,
        :misc
    ]
    specs = NamedTuple()
    for (key, value) in zip(keys(spec_dict), values(spec_dict))
        if key in spec_types
            spec_func = eval(Symbol("generate_$(key)_specs"))
            specs = merge(specs, [key => spec_func(value)])
        else 
            @error "Specified top-level component $key not part of HS3!"
        end
    end
    return specs
end

function generate_specs(path_to_file::String)
    generate_specs(file_to_dict(path_to_file))
end


"""
generate_misc_specs(dict)

Generate miscellaneous dictionaries from the provided dictionary.

# Arguments
- `dict::AbstractDict`: The dictionary containing the miscellaneous specifications.

# Returns
The provided dictionary as-is.

"""
generate_misc_specs(dict::AbstractDict) = Dict(dict)