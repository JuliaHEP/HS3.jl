# This file is a part of HS3.jl, licensed under the MIT License (MIT).
"""
    struct HistfactPDF

The `HistfactPDF` struct represents a probability density function (PDF) defined within the context of a HistFactory model given by LiteHF.jl.

# Fields
- `channel::Tuple`: A tuple representing the channel information associated with the PDF.
- `prior::NamedTuple`: A named tuple representing the prior/modifier information associated with the PDF.

"""
struct HistfactPDF 
    channel::Function 
    prior::NamedTuple 
    order::Vector
    customs::NamedTuple
end






"""
'_create_custom_nt' 

Support function for 'custom' modifiers 
"""
function _create_custom_nt(customs::AbstractArray)
    if length(customs) > 0
        list = [(first(item)) for item in customs]
        return NamedTupleTools.namedtuple(list)
    else
        return (;)
    end
end


"""
    make_histfact(histfact_spec::HistFactorySpec,  name::Symbol, function_specs=(;))

Constructs a HistfactPDF object based on the provided HistFactorySpec and channel name using the LiteHF.jl package. 
The specifiactions are first restructured to accomodate for the structure needed by LiteHF.jl which is then used to build a HistFactory object.

# Arguments
- `histfact_spec::HistFactorySpec`: The HistFactorySpec containing the specification for the HistFactory model.
- `name::Symbol`: The name of the channel associated with the HistFactoryPDF.
-`function_specs`: NamedTuple for function_specs/distribution_specs; association with functions is neccesary for custom modifier and constraint terms in normfactors

# Returns
A HistfactPDF object representing the constructed probability density function.

"""

function make_histfact(histfact_spec::HistFactorySpec, name::Symbol, function_specs=(;))
    channel_dict  = build_channel(histfact_spec, name, function_specs)
    channel = _create_expected(channel_dict) # the "pdf" representing the expected yields 
    custom = _create_custom_nt(channel_dict[4]) # an object representing the custom modifier
    order = setdiff(channel[2], collect(keys(custom))) # match parameters to the correct modifiers
    HistfactPDF(channel[1], channel[3], order, custom)
end

