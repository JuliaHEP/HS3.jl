# This file is a part of HS3.jl, licensed under the MIT License (MIT).

import Test
import Aqua


Test.@testset "Aqua tests" begin
    Aqua.test_all(HS3)
end # testset
