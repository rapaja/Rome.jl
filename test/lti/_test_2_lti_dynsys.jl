@testset "LTI dynsys tests" begin

    N = 10
    types = [Int32, Int64, Float32, Float64, ComplexF32, ComplexF64]
    atol = 1e-2
    
    @testset "differentiator tests" begin
        Δt = 0.1
        s = Romeo.LTI.Diff{Float64}(1)

        k = Romeo.LTI.kernel(s, N, Δt)
        @test k[1] == 1/Δt || "The first element of the differentiation kernel should be equal to 1/Δt"
        @test k[2] == -1/Δt || "The second element of the differentiation kernel should be -1/Δt"
        for i = 3:N
            @test k[i] == 0 || "All elements with index ≥3 of the differentiation kernel should be zero."
        end # for

        ys = Romeo.LTI.step_resp(s, N, Δt)
        @test ys[1] == 1/Δt || "The first element of the step response of a differentiator should be 1/Δt."
        for i = 2:N
            @test ys[i] == 0 || "All elements of the step response of a differentiator should be zero, except the first one."
        end # for        
    end # testset

    @testset "identity operator tests" begin
        Δt = 0.1
        s = Romeo.LTI.Diff{Float64}(0)

        k = Romeo.LTI.kernel(s, N, Δt)
        @test k[1] == 1 || "The first element of the identity kernel should be equal to 1"
        for i = 2:N
            @test k[i] == 0 || "All elements with index ≥2 of the identity kernel should be zero."
        end # for   

        ys = Romeo.LTI.step_resp(s, N, Δt)
        for i = 1:N
            @test ys[i] == 1 || "All elements of the step response of the identity kernel should be equal to 1."
        end # for 
    end # testset

    @testset "integrator tests" begin
        Δt = 1
        sinv = Romeo.LTI.Diff{Float64}(-1)
        k = Romeo.LTI.kernel(sinv, N, Δt)
        for i = 1:N
            @test k[i] == Δt || "All elements of the integration kernel should be equal to Δt."
        end # for

        ys = Romeo.LTI.step_resp(sinv, N, Δt)
        for i = 1:N
            @test ys[i] == i*Δt || "All elements of the integration kernel should be equal to i*Δt."
        end # for
    end # testset

    @testset "fractional integrator tests" begin
        Δt = 1
        int_0_5 = Romeo.LTI.Diff{Float64}(-0.5)
        int_1_5 = Romeo.LTI.Diff{Float64}(-0.5)

        ys = Romeo.LTI.impulse_resp(int_0_5, N, Δt)
        yi = Romeo.LTI.impulse_resp(int_1_5, N, Δt)
        
        @test all(ys .== yi) || "Step response of half integrator should match impulse response of 1.5 integrator."
    end

end # testset