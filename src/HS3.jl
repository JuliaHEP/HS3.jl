# This file is a part of HS3.jl, licensed under the MIT License (MIT).

"""
    HS3

Template for Julia packages.
"""

module HS3
import JSON3, PropDicts, Distributions, ValueShapes, StatsBase, DensityInterface, LiteHF, Parameters
import NamedTupleTools
abstract type AbstractHS3Spec end

include("BuiltinDist.jl")
include("BuiltinAxes.jl")
include("BuiltinDomains.jl")
include("BuiltinParameterPoints.jl")
include("read_in.jl")
include("BuiltinFunctions.jl")
include("MetaInformation.jl")
include("BuiltinData.jl")
include("TopologicalSort.jl")
include("BuiltinHistFact.jl")
include("BuiltinLikelihood.jl")
include("BuiltinAnalyses.jl")
include("Validator.jl")

#include("BuiltinOutput.jl")
end # module
