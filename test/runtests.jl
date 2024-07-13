using OrbitPropagationLibrarySOFA
using Test
using StaticArrays

const tests = [
    "Timing/JDate",
    "Timing/Conversions",
    "Timing/gmst",
    "Coordinates/rotations",
]
@testset "OrbitPropagationLibrarySOFA.jl" begin
    @testset "Test $t" for t in tests
        include("$t.jl")
    end
end
