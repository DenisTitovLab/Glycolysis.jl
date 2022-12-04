using Glycolysis
using Test

@testset "Glycolysis.jl" begin
    # Write your tests here.
    @test eltype(glycolysis_params) == Float64
    @test eltype(glycolysis_init_conc) == Float64
end
