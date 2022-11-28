# This file is a part of HS3.jl, licensed under the MIT License (MIT).

using Test
using HS3
import Documenter

Documenter.DocMeta.setdocmeta!(
    HS3,
    :DocTestSetup,
    :(using HS3);
    recursive=true,
)
Documenter.doctest(HS3)
