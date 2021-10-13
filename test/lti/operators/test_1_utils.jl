@testset "operators.utils tests" begin

    N = 10
    types = [Int32, Int64, Float32, Float64, ComplexF32, ComplexF64] 
    
    @testset "convolution/deconvolution tests" begin
        for T in types
            if T <: Integer
                g = rand(-10 * one(T):10 * one(T), N)
                h = rand(-10 * one(T):10 * one(T), N)
            else
                g = rand(T, N)
                h = rand(T, N)
            end
            f = Romeo.LTI.Operators.convolve(g, h)
            ga = Romeo.LTI.Operators.deconvolve(f, h)
            @test all(g .≈ ga) || "Error for type = $T"
        end # for
    end # testset

    @testset "convroot tests / positive first element" begin
        for T in types
            if T <: Integer
                g = rand(-10 * one(T):10 * one(T), N)
                h = Romeo.LTI.Operators.convolve(g, g)
            else
                g = rand(T, N)
                h = Romeo.LTI.Operators.convolve(g, g)
            end # if
            ga = Romeo.LTI.Operators.convroot(h)
            @test all(g .≈ ga) || "Error for type = $T"
        end # for        
    end # testset

    @testset "convroot tests / negative first element" begin
        for T in types
            if T <: Integer
                h = rand(-10 * one(T):10 * one(T), N)
                h[1] = rand(1 * one(T):10 * one(T))
            else
                h = rand(T, N)
                h[1] = abs(h[1]) + 0.01
            end # if
            g = Romeo.LTI.Operators.convroot(h)
            ha = Romeo.LTI.Operators.convolve(g, g)
            @test all(h .≈ ha) || "Error for type = $T"
        end # for        
    end # testset

end # testset