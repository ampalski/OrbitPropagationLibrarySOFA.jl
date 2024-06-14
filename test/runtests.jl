using OrbitPropagationLibrarySOFA
using Test

const tests = [
    "Timing/JDate",
    "Timing/Conversions",
]
@testset "OrbitPropagationLibrarySOFA.jl" begin
    @testset "Test $t" for t in tests
        include("$t.jl")
    end
end
