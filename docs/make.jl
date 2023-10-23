# Use
#
#     DOCUMENTER_DEBUG=true julia --color=yes make.jl local [nonstrict] [fixdoctests]
#
# for local builds.

using Documenter
using MJDSigGen

# Doctest setup
DocMeta.setdocmeta!(
    MJDSigGen,
    :DocTestSetup,
    :(using MJDSigGen);
    recursive=true,
)

makedocs(
    sitename = "MJDSigGen",
    modules = [MJDSigGen],
    format = Documenter.HTML(
        prettyurls = !("local" in ARGS),
        canonical = "https://legend-exp.github.io/MJDSigGen.jl/stable/"
    ),
    pages = [
        "Home" => "index.md",
        "API" => "api.md",
        "LICENSE" => "LICENSE.md",
    ],
    doctest = ("fixdoctests" in ARGS) ? :fix : true,
    linkcheck = !("nonstrict" in ARGS),
    warnonly = ("nonstrict" in ARGS),
)

deploydocs(
    repo = "github.com/legend-exp/MJDSigGen.jl.git",
    forcepush = true,
    push_preview = true,
)
