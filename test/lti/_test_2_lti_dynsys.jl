@testset "LTI dynsys tests" begin

    N = 10
    types = [Int32, Int64, Float32, Float64, ComplexF32, ComplexF64]
    atol = 1e-2
    
    @testset "Diff tests" begin
        s = Romeo.LTI.Diff{Float64}(1)
        k = Romeo.LTI.kernel(s, N)
        @test k[1] == 1 && k[2] == -1 || "Wrong kernel."
    end # testset

end # testset