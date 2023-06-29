# This file is a part of HS3.jl, licensed under the MIT License (MIT).

"""
    HS3

Template for Julia packages.
"""

module HS3
    import JSON3, Distributions, ValueShapes, StatsBase, DensityInterface, LiteHF, NamedTupleTools
    using Parameters 

    abstract type AbstractHS3Spec end

    include("read_in.jl")
    include("utils.jl")
   
    include("specs/otherspecs.jl")
    include("specs/parameter_pointspecs.jl")
    include("specs/spec_utils.jl")
    include("specs/domainspecs.jl")
    include("specs/distributionspecs.jl")
    include("specs/functionspecs.jl")
    include("specs/dataspecs.jl")
    include("specs/histfactspecs.jl")
    include("specs/likelihoodspecs.jl")
    include("specs/analysesspecs.jl")  

    for file in readdir("src/builders", sort = false)
        include(joinpath("builders/", file))
    end
    include("Validator.jl")
    #include("output/specs_to_dict.jl")
end 
