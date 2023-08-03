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
const spec_types = [
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

function _gen_specs(spec_dict::Vector)
    specs = (;)
    names = [element[1] for element in spec_dict]
    specs = [eval(Symbol("generate_$(element[1])_specs"))(element[2]) for element in spec_dict]
    @assert all(x -> x in spec_types, names)
    #for element in spec_dict
    #    key = element[1]
    #    value = element[2]
#
    #    if key in spec_types
    #        @info "The specs of field $(key) are being made"
    #        spec_func = eval(Symbol("generate_$(key)_specs"))
    #        specs = merge(specs, [key => spec_func(value)])
    #        @info "The specs of field $(key) have been made"
    #    else 
    #        @error "Specified top-level component $key not part of HS3!"
    #    end
    #end
    return (; zip(names, specs)...)
end

function generate_specs(spec_dict::AbstractDict)
    chunks = Iterators.partition(spec_dict, length(keys(spec_dict)) ÷ Threads.nthreads()+1)
    tasks = map(chunks) do chunk
        Threads.@spawn _gen_specs(chunk)
    end
    chunk_sums = fetch.(tasks)
    return merge(chunk_sums...)
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