# Use
#
#     DOCUMENTER_DEBUG=true julia --color=yes make.jl local [nonstrict] [fixdoctests]
#
# for local builds.

using Documenter
using HS3

# Doctest setup
DocMeta.setdocmeta!(
    HS3,
    :DocTestSetup,
    :(using HS3);
    recursive=true,
)

makedocs(
    sitename = "HS3",
    modules = [HS3],
    format = Documenter.HTML(
        prettyurls = !("local" in ARGS),
        canonical = "https://oschulz.github.io/HS3.jl/stable/"
    ),
    pages = [
        "Home" => "index.md",
        "API" => "api.md",
        "LICENSE" => "LICENSE.md",
    ],
    doctest = ("fixdoctests" in ARGS) ? :fix : true,
    linkcheck = !("nonstrict" in ARGS),
    strict = !("nonstrict" in ARGS),
)

deploydocs(
    repo = "github.com/oschulz/HS3.jl.git",
    forcepush = true,
    push_preview = true,
)
