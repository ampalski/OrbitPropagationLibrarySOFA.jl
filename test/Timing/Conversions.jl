@testset "Time Conversions Vallado" begin
    j, m = dateVec2JDate([2004.0, 5, 14, 16, 43, 0.0], isUTC=true)
    m_UT1 = UTC2UT1(m, type=:MJD)
    m_TAI = UTC2TAI(m, type=:MJD)
    m_TT = UT12TT(m_UT1, type=:MJD)
    m_TDB = TT2TDB(m_TT, type=:MJD)

    d = JDate2dateVec(m_UT1, type=:MJD)
    @test d[5] == 42
    @test isapprox(d[6], 59.5367, atol=0.0001)

    d = JDate2dateVec(m_TAI, type=:MJD)
    @test d[5] == 43
    @test isapprox(d[6], 32.0, atol=0.0001)

    d = JDate2dateVec(m_TT, type=:MJD)
    @test d[5] == 44
    @test isapprox(d[6], 4.184, atol=0.0001)

    d = JDate2dateVec(m_TDB, type=:MJD)
    @test d[5] == 44
    @test isapprox(d[6], 4.1856, atol=0.001)
    # Slightly less tolerance due to the approximation made in TDB conversions
end

# From t_sofa_c.c
@testset "Time Conversions SOFA" begin
    #UTC2TAI
    j = [2453750.5, 0.892100694]
    jtai = UTC2TAI(j, type=:JD)
    @test jtai[1] == 2453750.5
    @test isapprox(jtai[2], 0.8924826384444444444, atol=1e-12)

    #UTC2UT1
    jut1 = UTC2UT1(j, type=:JD)
    @test jut1[1] == 2453750.5
    @test isapprox(jut1[2], 0.8921045609398149, atol=1e-12)
    # Slightly different than t_sofa_c to account for more precise Δut1 term

end
