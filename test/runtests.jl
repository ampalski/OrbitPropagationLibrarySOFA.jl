using OrbitPropagationLibrarySOFA
using Test
using StaticArrays
using Aqua
using JET

const tests = [
    "Timing/JDate",
    "Timing/Conversions",
    "Timing/gmst",
    "Coordinates/rotations",
    "Extras/aqua",
    "Extras/jet",
]
@testset "OrbitPropagationLibrarySOFA.jl" begin
    @testset "Test $t" for t in tests
        include("$t.jl")
    end
end
