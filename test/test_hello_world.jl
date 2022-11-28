# This file is a part of HS3.jl, licensed under the MIT License (MIT).

using HS3
using Test


@testset "hello_world" begin
    @test HS3.hello_world() == 42
end
