@testset "LTI system kernels tests" begin
    N = 10
    types = [Int32, Int64, Float32, Float64, ComplexF32, ComplexF64]
    Δt = 0.1
    atol = Δt

    @testset "Diff tests" begin

        @testset "identity operator Diff(0) tests" begin
            id = Romeo.LTI.Diff(0)
            k = Romeo.LTI.kernel(id, N, Δt)
            @test k[1] == 1 || "The first element of the identity kernel should be equal to 1"
            for i = 2:N
                @test k[i] == 0 || "All elements with index ≥2 of the identity kernel should be zero."
            end # for   

            ys = Romeo.LTI.step_resp(id, N, Δt)
            for i = 1:N
                @test ys[i] == 1 || "All elements of the step response of the identity operator should be equal to 1."
            end # for 
        end # testset

        @testset "first order integral operator tests" begin
            int = Romeo.LTI.Diff(-1)
            k = Romeo.LTI.kernel(int, N, Δt)
            for i = 1:N
                @test k[i] == Δt || "All elements of the first order integral kernel should be equal to Δt."
            end # for  

            yi = Romeo.LTI.impulse_resp(int, N, Δt)
            for i = 1:N
                @test yi[i] == 1 || "All elements of the impulse response of the first order integrator should be equal to 1."
            end # for

            ys = Romeo.LTI.step_resp(int, N, Δt)
            for i = 1:N
                expected = (i - 1) * Δt
                @test isapprox(ys[i], expected; atol=atol) || "Step response of the first order integrator at index $i : $(ys[i]) ≠ $expected"
            end # for  
        end # testset

        @testset "second order integral operator tests" begin

            int2 = Romeo.LTI.Diff(-2)
            k = Romeo.LTI.kernel(int2, N, Δt)
            for i = 1:N
                expected = 1 // 2 * Δt^2
                @test isapprox(k[i], expected; atol=atol) || "Element of the second order integral kernel at index $i : $(k[i]) ≠ $expected."
            end # for  

        end # testset

        @testset "differentiator tests" begin
            Δt = 0.1
            s = Romeo.LTI.Diff{Float64}(1)

            k = Romeo.LTI.kernel(s, N, Δt)
            @test k[1] == 1 / Δt || "The first element of the differentiation kernel should be equal to 1/Δt"
            @test k[2] == -1 / Δt || "The second element of the differentiation kernel should be -1/Δt"
            for i = 3:N
                @test k[i] == 0 || "All elements with index ≥3 of the differentiation kernel should be zero."
            end # for

            ys = Romeo.LTI.step_resp(s, N, Δt)
            @test ys[1] == 1 / Δt || "The first element of the step response of a differentiator should be 1/Δt."
            for i = 2:N
                @test ys[i] == 0 || "All elements of the step response of a differentiator should be zero, except the first one."
            end # for        
        end # testset

    end

end # testset