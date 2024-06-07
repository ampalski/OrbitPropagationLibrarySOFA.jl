@testset "JDate" begin
    date1 = [2024, 5, 29, 8.0, 40, 30]
    date2 = [1999.0, 1, 1, 0, 0, 0]

    JD, MJD = dateVec2JDate(date1)
    JD2, MJD2 = dateVec2JDate(date2)

    @test isapprox(sum(JD), 2460459.8614598; atol=1e-5)
    @test isapprox(sum(JD2[1]), 2451179.5000000; atol=1e-5)
    @test JD2[2] == 0.0
    @test MJD2[2] == 0.0
    @test MJD2[1] == 51179.0

end
