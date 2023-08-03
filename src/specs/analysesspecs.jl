# This file is a part of HS3.jl, licensed under the MIT License (MIT).

"""
    @with_kw struct AnalysesSpec <: AbstractHS3Spec
        likelihood::String
        parameter_domain::String
        parameters_of_interest::Union{AbstractArray{String}, Nothing} = nothing
        parameter_init::Union{String, Nothing} = nothing
        prior::Union{String, Nothing} = nothing
    end
Struct representing specifications for analyses in HS3.

Fields:
- `likelihood`: String specifying the likelihood.
- `parameter_domain`: String specifying the domain.
- `parameters_of_interest`: Optional array of Strings specifying the parameters of interest.
- `parameter_init`: Optional String specifying the parameter initialization.
- `prior`: Optional String specifying the prior distribution.
"""
@with_kw struct AnalysesSpec <: AbstractHS3Spec
    likelihood::String 
    domains::Union{AbstractArray{String}, Nothing, AbstractArray, String} = nothing
    parameters_of_interest::Union{AbstractArray{String}, Nothing} = nothing
    parameter_init::Union{String, Nothing} = nothing
    prior::Union{String, Nothing} = nothing
end


"""
    generate_analyses_specs(analyses_array::AbstractArray)

Generate analysis specifications from an array of analysis data.

This function takes an array of analysis data and generates analysis specifications using the `AnalysesSpec` struct. 
Each element in the `analyses_array` represents an analysis, and the function creates a named tuple of analysis specifications
where the analysis name is used as the field name.

Arguments:
- `analyses_array`: An abstract array containing analysis data.

Returns:
- NamedTuple: A named tuple of analysis specifications.

Example:
```julia
analyses_array = [...]  # Array of analysis data
specs = generate_analyses_specs(analyses_array)
"""
function generate_analyses_specs(analyses_array::AbstractArray)
    #specs = NamedTuple{}()
    names = [Symbol(analysis.name) for analysis in analyses_array]
    specs = [AnalysesSpec(; NamedTuple(filter(entry -> entry[1] != :name, analysis))...) for analysis in analyses_array]
    #for analysis in analyses_array
    #    name = Symbol(analysis.name)
    #    data = NamedTuple(filter(entry -> entry[1] != :name, analysis))
    #    specs = merge(specs, (name => AnalysesSpec(; data...),),)
    #end
    #specs
    return (; zip(names, specs)...)
end
