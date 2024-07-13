@testset "GMST Vallado" begin
    date1 = [2004.0, 4, 6, 7, 51, 27.946047]

    JD, MJD = datevec2jdate(date1, system=:UT1)

    @test isapprox(gmst(JD; model=82), 5.459562588350443; atol=1e-8)
    @test isapprox(gast(JD; model=94), 5.459507978743478; atol=1e-8)

end

@testset "GMST SOFA" begin
    j = JDate(SA[2400000.5, 53736.0], :UT1)

    @test isapprox(gmst(j; model=82), 1.754174981860675096; atol=1e-12)
    @test isapprox(gast(j; model=94), 1.754166136020645203; atol=1e-9)
    # 1e-9 vs 1e-12 difference for GAST is as simple as using TT for eqeq
    # instead of UT1 ... spec calls for eqeq being calculated with TT, but the
    # SOFA test uses UT1
end

