@testset "LTI utils tests" begin

    N = 10
    types = [Int32, Int64, Float32, Float64, ComplexF32, ComplexF64]
    atol = 1e-2
    
    @testset "convolution/deconvolution tests" begin
        zerotol = 5
        for T in types
            if T <: Integer
                g = rand(-10 * one(T):10 * one(T), N)
                h = rand(-10 * one(T):10 * one(T), N)
                if h[1] == 0
                    h[1] = convert(T, rand([-1, 1]))
                end
            else
                g = rand(T, N)
                h = rand(T, N)
                if abs(h[1]) < zerotol
                    if T <: Real
                        h[1] = convert(T, zerotol) * sign(h[1])
                    else
                        h[1] = convert(T, zerotol * exp(angle(h[1]) * im))
                    end
                end
            end
            f = Romeo.LTI.convolve(g, h)
            ga = Romeo.LTI.deconvolve(f, h)
            @test isapprox(g, ga; atol=atol) || "Error for type = $T. max abs.err. = $(max(abs.(g - ga)...))"
        end # for
    end # testset

    @testset "convroot tests / positive first element" begin
        for T in types
            if T <: Integer
                g = rand(-10 * one(T):10 * one(T), N)
                g[1] = rand(1 * one(T):10 * one(T))
                h = Romeo.LTI.convolve(g, g)
            else
                g = rand(T, N)
                g[1] = abs(g[1]) + 0.01
                h = Romeo.LTI.convolve(g, g)
            end # if
            ga = Romeo.LTI.convroot(h)
            @test isapprox(g, ga; atol=atol) || "Error for type = $T. max abs.err. = $(max(abs.(g - ga)...))"
        end # for        
    end # testset

    @testset "convroot tests / negative first element" begin
        for T in types
            if T <: Integer
                h = rand(-10 * one(T):10 * one(T), N)
                h[1] = -rand(1 * one(T):10 * one(T))
            else
                h = rand(T, N)
                h[1] = -(abs(h[1]) + 0.01)
            end # if
            g = Romeo.LTI.convroot(h)
            ha = Romeo.LTI.convolve(g, g)
            @test isapprox(h, ha; atol=atol) || "Error for type = $T. max abs.err. = $(max(abs.(h - ha)...))"
        end # for        
    end # testset

end # testset