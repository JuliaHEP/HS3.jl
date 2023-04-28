# This file is a part of HS3.jl, licensed under the MIT License (MIT).

import Test


Test.@testset "HS3" begin
    #include("test_aqua.jl")
    #include("test_hello_world.jl")
    include("test_docs.jl")
end # testset
using HS3
HS3.hello_world()