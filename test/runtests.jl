using BiophysicalParameters
using Test
using SafeTestsets
using Aqua

@testset "BiophysicalParameters.jl" begin
    Aqua.test_all(BiophysicalParameters; ambiguities=false)
    @safetestset "provenance" include("provenance.jl")
    @safetestset "radiation"  include("radiation.jl")
end
