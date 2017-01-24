# This file is a part of MJDSigGen, licensed under the MIT License (MIT).

module MJDSigGen

__precompile__(false)


const deps_jl = joinpath(dirname(@__FILE__), "..", "deps", "deps.jl")
if isfile(deps_jl)
    include(deps_jl)
else
    error("MJDSigGen is not properly installed. Run Pkg.build(\"MJDSigGen\").")
end


include("types.jl")


end # module
