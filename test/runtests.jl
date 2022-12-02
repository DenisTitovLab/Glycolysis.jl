using Glycolysis
using Test

@testset "Glycolysis.jl" begin
    # Write your tests here.
    @test eltype(model_params) == Float64
    @test eltype(initial_concentrations) == Float64
end
