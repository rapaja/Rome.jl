@testset "LTI systems tests" begin
        
    s = Romeo.LTI.Diff(1)
    G = (s+1)/(s+2)
    zerosys = Romeo.LTI.ZeroSys()
    idsys = Romeo.LTI.UnitSys()

    @testset "Diff tests" begin
        @test s.α == 1
        @test Romeo.LTI.Diff(0) == idsys
    end

    @testset "Unary operators tests" begin
        @test +s == s
        @test -s == Romeo.LTI.ScaledSystem(-1, s)
        @test -(-s) == s
    end

    @testset "Series connection tests" begin
        @test 1*s == s*1 == s
        @test 0*s == s*0 == zerosys
        @test idsys * s == s * idsys == s
        @test zerosys * s == s * zerosys == zerosys
        @test 2*s == s*2 == Romeo.LTI.ScaledSystem(2, s)
        @test 4*(s*2) == Romeo.LTI.ScaledSystem(8, s)
        @test -(8*s) == Romeo.LTI.ScaledSystem(-8, s)
        @test (2*s)*0.5 == s
        @test (2*s)*G == s*(2*G) == 2*(s*G)
        @test (2*s)*(4*s) == 8*(s*s)
    end

    @testset "Parallel connection tests" begin
        @test s + s == 2*s
        @test 2*s + 4*s == 6*s
        @test s + 2*s == 2*s + s == 3*s
        @test s - s == zerosys
        @test s + zerosys == zerosys + s == s
        @test s + 0 == 0 + s == s
    end


    @testset "SquareRoot- and Power- systems tests" begin
        
    end #testset

end # testset
