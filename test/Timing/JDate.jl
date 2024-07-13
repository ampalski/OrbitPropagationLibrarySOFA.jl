@testset "JDate" begin
    date1 = [2024, 5, 29, 8.0, 40, 30]
    date2 = [1999.0, 1, 1, 0, 0, 0]

    JD, MJD = dateVec2JDate(date1)
    JD2, MJD2 = dateVec2JDate(date2)

    @test isapprox(sum(JD.epoch), 2460459.8614598; atol=1e-5)
    @test isapprox(sum(JD2.epoch[1]), 2451179.5000000; atol=1e-5)
    @test JD2.epoch[2] == 0.0
    @test MJD2.epoch[2] == 0.0
    @test MJD2.epoch[1] == 51179.0

    date1_new = JDate2dateVec(JD)
    date2_new = JDate2dateVec(JD2)

    @test date2_new == [1999.0, 1, 1, 0, 0, 0]
    @test date1_new[3] == 29.0
    @test isapprox(date1_new[6], 30.0; atol=1e-10)

end

@testset "JDATE SOFA" begin
    date1 = [1994.0, 6, 30, 23, 59, 60.13599]
    date2 = [2003.0, 6, 1, 0, 0, 0]
    JD, MJD = dateVec2JDate(date1, system=:UTC)
    JD2, MJD2 = dateVec2JDate(date2)

    @test isapprox(sum(JD.epoch), 2449534.49999, atol=1e-6)
    @test MJD2.epoch[1] == 52791.0
    @test MJD2.epoch[2] == 0

    date1_new = JDate2dateVec(JD)
    date2_new = JDate2dateVec(JD2)

    @test date2_new == [2003.0, 6, 1, 0, 0, 0]
    @test date1_new[3] == 30.0
    @test isapprox(date1_new[6], 60.13599; atol=1e-4)
end

@testset "fixJDate" begin
    date1 = [2024.0, 5, 29, 25, 40, 30]
    date2 = [2024.0, 5, 31, 24, 40, 30]
    date3 = [2024.0, 12, 31, 23, 59, 5000]
    date4 = [2023.0, 5, 31, 24, 40, 30]
    date5 = [2024.0, 6, -1, 12, 12, 12]
    date6 = [2024.0, 6, 1, -12, 12, 12]

    fixDateVec!(date1)
    fixDateVec!(date2)
    fixDateVec!(date3)
    fixDateVec!(date4)
    fixDateVec!(date5)
    fixDateVec!(date6)

    @test date1 == [2024.0, 5, 30, 1, 40, 30]
    @test date2 == [2024.0, 6, 1, 0, 40, 30]
    @test date3 == [2025.0, 1, 1, 1, 22, 20]
    @test date4 == [2023.0, 6, 1, 0, 40, 30]
    @test date5 == [2024.0, 5, 30, 12, 12, 12]
    @test date6 == [2024.0, 5, 31, 12, 12, 12]

end



