@testset "LTI systems tests" begin

    s = Romeo.LTI.Diff(1)
    H = 1 / (s + 1)
    G = (s + 1) / (s + 2)

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
        for sys in [s, G]
            @test 1 * sys == sys * 1 == sys
            @test 1.0 * sys == sys * 1.0 == sys
            @test 2 * sys == 2.0 * sys == sys * 2.0 == sys * 2
            @test 0 * sys == sys * 0 == 0.0 * sys == sys * 0.0 == zerosys
            @test -0 * sys == -sys * 0 == -0.0 * sys == -sys * 0.0 == zerosys
            @test idsys * sys == sys * idsys == sys
            @test zerosys * sys == sys * zerosys == zerosys
            @test 2 * sys == sys * 2 == Romeo.LTI.ScaledSystem(2, sys)
            @test 4 * (sys * 2) == Romeo.LTI.ScaledSystem(8, sys)
            @test -(8 * sys) == Romeo.LTI.ScaledSystem(-8, sys)
            @test (2 * sys) * 0.5 == sys
            @test (2 * sys) * (4 * sys) == 8 * (sys * sys)
        end
        @test (2 * s) * G == s * (2 * G) == 2 * (s * G)
    end

    @testset "Parallel connection tests" begin
        for sys in [s, G, H]
            @test sys + sys == 2 * sys
            @test 2 * sys + 4 * sys == 6 * sys
            @test sys + 2 * sys == 2 * sys + sys == 3 * sys
            @test sys - sys == zerosys
            @test sys + zerosys == zerosys + sys == sys
            @test sys + 0 == 0 + sys == sys + 0.0 == 0.0 + sys == sys
        end
    end

    @testset "Power system tests" begin
        orders = [0.5, 2, 2.5, 3]
        @test idsys^0 == idsys^0.0 == idsys
        @test zerosys^0 == zerosys^0.0 == idsys
        @test sqrt(idsys) == idsys
        @test sqrt(zerosys) == zerosys
        @test sqrt(s^0.5) == s^0.25
        for α in orders
            @test idsys^α == idsys
            @test zerosys^α == zerosys
        end
        for sys in [s, G]
            @test sys^0 == sys^0.0 == idsys
            @test sys^1 == sys^1.0 == sys
            @test sys^0.5 == sqrt(sys)
            @test sqrt(4 * sys) == 2 * sqrt(sys)
            for α in orders
                @test (2 * sys)^α == 2^α * sys^α
            end
            @test (sys^2)^3 == sys^6
        end
    end

    @testset "Rational systems tests" begin
        for sys in [s, G, H]
            @test sys / 1 == sys / 1.0 == sys == 1 \ sys == 1.0 \ sys
            @test sys / idsys == sys == idsys \ sys
            # @test_throws ErrorException sys / zerosys
            # @test_throws ErrorException zerosys \ sys
            @test sys / sys == idsys
            @test sys^(15 // 100) / sys^(1 // 10) == sys^(5 // 100)
            @test sys / sys^0.4 == sys^0.6
            @test sys^(12 // 10) / sys == sys^(2 // 10)
        end
        @test s / (4 * G) == 1 / 4 * s / G
        @test (6 * s) / G == 6 * (s / G)
    end


end # testset
