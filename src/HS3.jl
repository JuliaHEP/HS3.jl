# This file is a part of HS3.jl, licensed under the MIT License (MIT).

"""
    HS3

WIP package for HS3 - HEP Statistics Serialization Standard
    Tests and Docs might be added later.
"""

module HS3
    import JSON3, Distributions, ValueShapes, StatsBase, DensityInterface, NamedTupleTools
    using Unrolled
    using SpecialFunctions: logabsgamma, loggamma, gamma, erf
    using LogExpFunctions: xlogy
    using Parameters 
    import Distributions: logpdf, pdf
    using QuadGK
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

    include("builders/build_analyses.jl")
    include("builders/modifier_utils.jl")
    include("builders/build_data.jl")
    include("builders/build_distributions.jl")
    include("builders/build_functionals.jl")
    include("builders/build_functions.jl")
    include("builders/build_histfactory.jl")
    include("builders/build_likelihoods.jl")
    include("builders/build_other.jl")
    include("builders/histfactory_utils.jl")
    include("builders/topological_sort.jl")
end 
