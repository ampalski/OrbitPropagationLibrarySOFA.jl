@testset "Time Conversions Vallado" begin
    j, m = datevec2jdate([2004.0, 5, 14, 16, 43, 0.0])
    m_UT1 = convert_jd(m, :UT1)
    m_TAI = convert_jd(m_UT1, :TAI)
    m_TT = convert_jd(m_UT1, :TT)
    m_TDB = convert_jd(m_UT1, :TDB)

    d = jdate2datevec(m_UT1)
    @test d[5] == 42
    @test isapprox(d[6], 59.5367, atol=0.0001)

    d = jdate2datevec(m_TAI)
    @test d[5] == 43
    @test isapprox(d[6], 32.0, atol=0.0001)

    d = jdate2datevec(m_TT)
    @test d[5] == 44
    @test isapprox(d[6], 4.184, atol=0.0001)

    d = jdate2datevec(m_TDB)
    @test d[5] == 44
    @test isapprox(d[6], 4.1856, atol=0.001)
    # Slightly less tolerance due to the approximation made in TDB conversions
end

# From t_sofa_c.c
@testset "Time Conversions SOFA" begin
    #UTC2TAI
    j = JDate(SA[2453750.5, 0.892100694], :UTC)
    jtai = convert_jd(j, :TAI)
    # jtai = UTC2TAI(j, type=:JD)
    @test jtai.epoch[1] == 2453750.5
    @test isapprox(jtai.epoch[2], 0.8924826384444444444, atol=1e-12)

    #UTC2UT1
    jut1 = convert_jd(j, :UT1)
    # jut1 = UTC2UT1(j, type=:JD)
    @test jut1.epoch[1] == 2453750.5
    @test isapprox(jut1.epoch[2], 0.8921045609398149, atol=1e-12)
    # Slightly different than t_sofa_c to account for more precise Δut1 term

    #TT2TAI
    jtt = JDate(SA[2453750.5, 0.892482639], :TT)
    jtai = convert_jd(jtt, :TAI)
    # jtai = TT2TAI(jtt)
    @test jtai.epoch[1] == 2453750.5
    @test isapprox(jtai.epoch[2], 0.892110139, atol=1e-12)

    #UT12UTC
    jut1 = JDate(SA[2453750.5, 0.892104561], :UT1)
    jutc = convert_jd(jut1, :UTC)
    # jutc = UT12UTC(jut1, type=:JD)
    @test jutc.epoch[1] == 2453750.5
    @test isapprox(jutc.epoch[2], 0.8921006940601852, atol=1e-12)
    # Slightly different than t_sofa_c to account for more precise Δut1 term

    #TAI2UT1
    jtai = JDate(SA[2453750.5, 0.892482639], :TAI)
    jut1 = convert_jd(jtai, :UT1)
    # jut1 = TAI2UT1(jtai, type=:JD)
    @test jut1.epoch[1] == 2453750.5
    @test isapprox(jut1.epoch[2], 0.8921045614953704, atol=1e-12)
    # Slightly different than t_sofa_c to account for more precise Δut1 term

    #UT12TT
    jut1 = JDate(SA[2453750.5, 0.892104561], :UT1)
    jtt = convert_jd(jut1, :TT)
    # jtt = UT12TT(jut1, type=:JD)
    @test jtt.epoch[1] == 2453750.5
    @test isapprox(jtt.epoch[2], 0.8928551385046295, atol=1e-12)

end
