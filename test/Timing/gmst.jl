@testset "GMST Vallado" begin
    date1 = [2004.0, 4, 6, 7, 51, 27.946047]

    JD, MJD = dateVec2JDate(date1)

    @test isapprox(GMST(JD; model=82), 5.459562588350443; atol=1e-8)
    @test isapprox(GAST(JD; model=94), 5.459507978743478; atol=1e-8)

end

@testset "GMST SOFA" begin
    j = [2400000.5, 53736.0]

    @test isapprox(GMST(j; model=82), 1.754174981860675096; atol=1e-12)
    @test isapprox(GAST(j; model=94), 1.754166136020645203; atol=1e-9)
    #TODO: figure out if the 1e-9 vs 1e-12 difference for GAST is as simple as
    #using TT for eqeq instead of UT1
end

