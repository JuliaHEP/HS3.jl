# This file is a part of HS3.jl, licensed under the MIT License (MIT).

function validate_file(path::String)
    #check if JSON-File is valid
    @info "Now checking validity of JSON"
    dict = try
        file_to_dict(path)
    catch e
        @error "Invalid JSON File"
        throw(e)
    end

    #check for HS3 version number
    @info "Now checking HS3 version number"
    try 
        @assert haskey(dict.metadata, :hs3_version)
        try 
            @assert dict.metadata.hs3_version == 0.2
        catch e
            @error "HS3 version declaration is not compatible with this HS3-implementation"
        end
    catch e
        @error "HS3 version declaration is missing"
    end

    #check for correct top-level components
    toplevel = [:distributions, :functions, :data, :likelihoods, :domains, :parameter_points, :analyses, :metadata, :misc]
    for key in keys(dict)
        @assert key in toplevel "The component <$key> is not part of HS3" 
    end

    #test for correct type implementations
    @info "Now testing for correct type implementations"
    dist_specs = NamedTuple()
    if haskey(dict, :distributions)   
        @assert typeof(dict.distributions) <: JSON3.Array "Component <distributions> has to be of Array-Type" 
        dist_specs = merge(dist_specs, 
        try 
            make_distspecs(dict.distributions)
        catch e
            @error "Distribution is not defined correctly"
            rethrow(e)
        end
        )
        @assert length(dist_specs) == length(dict.distributions) "All distributions need to have a unique name"
    else
        @info "Component <distributions> not found"
    end
    func_specs = NamedTuple()
    if haskey(dict, :functions)   
        @assert typeof(dict.functions) <: JSON3.Array "Component <functions> has to be of Array-Type" 
        func_specs = merge(func_specs, 
        try 
            make_funcspecs(dict.functions)
        catch e
            @error "Function is not defined correctly"
            rethrow(e)
        end
        )
        @assert length(func_specs) == length(dict.functions) "All functions need to have an unique name"
    else
        @info "Component <function> not found"
    end

    data_specs = NamedTuple()
    if haskey(dict, :data)   
        @assert typeof(dict.data) <: JSON3.Array "Component <data> has to be of Array-Type" 
        data_specs = merge(data_specs, 
        try 
            make_dataspecs(dict.data)
        catch e
            @error "Data is not defined correctly"
            rethrow(e)
        end
        )
        @assert length(data_specs) == length(dict.data) "All data sets need to have an unique name"
    else
        @info "Component <data> not found"
    end

    likelihoods = NamedTuple()
    if haskey(dict, :likelihoods)
        @assert haskey(dict, :distributions) "Likelihoods depend on distributions, distributions not defined"
        @assert haskey(dict, :data) "Likelihoods depend on data, data not defined"
        @assert typeof(dict.likelihoods) <: JSON3.Array "Component <likelihoods> has to be of Array-Type" 
        for i in 1:length(dict.likelihoods)
            @assert haskey(dict.likelihoods[i], :name) "Likelihood $i is missing a <name>"
            @assert haskey(dict.likelihoods[i], :distributions) "Likelihood <$(dict.likelihoods[i].name)> is missing the key <distributions>"
            @assert haskey(dict.likelihoods[i], :data) "Likelihood <$(dict.likelihoods[i].name)> is missing the key <data>"
            @assert length(dict.likelihoods[i].distributions) == length(dict.likelihoods[i].data) "Likelihood <$(dict.likelihoods[i].name)> needs to have the same number of <distributions> and <data> entries"
        end
        likelihoods = make_likelihood_specs(dict.likelihoods)
    else 
        @info "Component <likelihoods> not found"
    end

    domain_specs = NamedTuple()
    if haskey(dict, :domains)   
        @assert typeof(dict.domains) <: JSON3.Array "Component <domains> has to be of Array-Type" 
        domain_specs = merge(domain_specs, 
        try 
            make_domainspecs(dict.domains)
        catch e
            @error "<Domains> is not defined correctly"
            rethrow(e)
        end
        )
        #println(domain_specs)
        #for i in 1:length(domain_specs)
        #    for j in 1:length(domain_specs[i])
        #        @assert typeof(domain_specs[i][j]) <: NamedTuple{(:max, :min), T} where T<:Tuple{<:Real, <:Real} "Parameter $(keys(domain_specs[i][j]) of domain $(keys(domain_specs)[i]) is not correctly defined "
        #    end
        #end
    else
        @info "Component <domains> not found"
    end
    point_specs = NamedTuple()
    if haskey(dict, :parameter_points)   
        @assert typeof(dict.parameter_points) <: JSON3.Array "Component <parameter_points> has to be of Array-Type" 
        point_specs = merge(point_specs, 
        try 
            make_points(dict.parameter_points)
        catch e
            @error "<parameter_points> is not defined correctly"
            rethrow(e)
        end
        )
    else
        @info "Component <parameter_points> not found"
    end
    @assert typeof(dict.metadata) <: AbstractDict "Component <metadata> has to be of Dictionary-Type" 
    check_metadata(dict.metadata )

    if haskey(dict, :misc)
        @assert typeof(dict.misc) <: AbstractDict  "Component <misc> has to be of Dictionary-Type" 
    else
        @info "Component <misc> not found"
    end

    if haskey(dict, :analyses)
        @assert typeof(dict.analyses) <: JSON3.Array "Component <analyses> has to be of Array-Type" 
        for ana in dict.analyses
            @assert typeof(ana) <: AbstractDict "Sub-components of <analyses> have to be of Dictionary-Type" 
            for comp in [:name, :likelihood, :parameter_domain, :data_domain]
                @assert haskey(ana, comp) "Top-level component <analyses> is missing sub-component <$(comp)>"
            end
            @assert (Symbol(ana.likelihood) in keys(likelihoods)) "Specified likelihood <$(ana.likelihood)> not found in analysis <$(ana.name)>"
            @assert (Symbol(ana.parameter_domain) in keys(domain_specs)) "Specified parameter domain <$(ana.parameter_domain)> not found in analysis <$(ana.name)>" 
            @assert (Symbol(ana.data_domain) in keys(domain_specs)) "Specified data domain <$(ana.parameter_domain)> not found in analysis <$(ana.name)>" 
            names = (; [:par_names => keys((domain_specs[Symbol(ana.parameter_domain)]).params), :var_names => keys((domain_specs[Symbol(ana.data_domain)]).params)]...)
            print(likelihoods[Symbol(ana.likelihood)])
            check_likelihood(likelihoods[Symbol(ana.likelihood)], dist_specs, names)
        end
    elseif haskey(dict, :likelihoods)
    
    
    end    
end

export validate_file

function check_likelihood(likelihood::HS3.LikelihoodSpec, dist_specs, names::NamedTuple)
    #likelihood.data 
    print(dist_specs[Symbol(likelihood.distributions)])
end

function check_metadata(dict::AbstractDict)
    allowed_keys = [:hs3_version, :packages, :authors, :publications, :description]
    if !haskey(dict, :hs3_version)
        @warn "HS3 version is missing. File/Content might be incompatible with this implementation"
    elseif dict.hs3_version != 0.2 
        @warn "Specified HS3 version $(dict.hs3_version) not compatible with HS3 version of this implementation (v. 0.2)"
    end
    for key in keys(dict)
        @assert (key in allowed_keys) "The component <$(key)> is not an allowed component in Metadata"
    end
end