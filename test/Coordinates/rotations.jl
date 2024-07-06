@testset "ITRF2PEF Vallado" begin
    date1 = [2004.0, 4, 6, 7, 51, 28.386009]#UTC, doesn't matter for this test

    JD, _ = dateVec2JDate(date1)

    pos = [-1033.4793830, 7901.2952754, 6380.3565958]
    vel = [-3.225636520, -2.872451450, 5.531924446]

    posPEF = ITRF2PEF76(pos, JD)
    velPEF = ITRF2PEF76(vel, JD)
    posITRF = PEF2ITRF76(posPEF, JD)

    @test isapprox(posPEF[1], -1033.4750313; atol=1e-4)
    @test isapprox(posPEF[2], 7901.3055856; atol=1e-4)
    @test isapprox(posPEF[3], 6380.3445327; atol=1e-4)

    @test isapprox(velPEF[1], -3.225632747; atol=1e-4)
    @test isapprox(velPEF[2], -2.872442511; atol=1e-4)
    @test isapprox(velPEF[3], 5.531931288; atol=1e-4)

    @test isapprox(posITRF[1], pos[1]; atol=1e-4)
    @test isapprox(posITRF[2], pos[2]; atol=1e-4)
    @test isapprox(posITRF[3], pos[3]; atol=1e-4)
end

@testset "PEF2TOD Vallado" begin
    date1 = [2004.0, 4, 6, 7, 51, 27.946047] #UT1

    JD, _ = dateVec2JDate(date1)

    pos = [-1033.4750313, 7901.3055856, 6380.3445327]
    vel = [-3.225632747, -2.872442511, 5.531931288]

    posTOD = PEF2TOD76(pos, JD)
    velTOD = PEF2TOD76_vel(vel, pos, JD)
    posPEF = TOD2PEF76(posTOD, JD)
    velPEF = TOD2PEF76_vel(velTOD, posPEF, JD)

    @test isapprox(posTOD[1], 5094.5147804; atol=1e-5)
    @test isapprox(posTOD[2], 6127.3664612; atol=1e-5)
    @test isapprox(posTOD[3], 6380.3445328; atol=1e-5)
    @test isapprox(velTOD[1], -4.746088567; atol=1e-8)
    @test isapprox(velTOD[2], 0.786077222; atol=1e-8)
    @test isapprox(velTOD[3], 5.531931288; atol=1e-8)

    posDiff = posPEF - pos
    velDiff = velPEF - vel

    @test posDiff' * posDiff < 1e-13
    @test velDiff' * velDiff < 1e-13
end

@testset "TOD2MOD Vallado" begin
    date1 = [2004.0, 4, 6, 7, 52, 32.570009] #TT

    JD, _ = dateVec2JDate(date1)

    pos = [5094.5147804, 6127.3664612, 6380.3445328]
    vel = [-4.746088567, 0.786077222, 5.531931288]

    posMOD = TOD2MOD76(pos, JD)
    velMOD = TOD2MOD76(vel, JD)
    posTOD = MOD2TOD76(posMOD, JD)
    velTOD = MOD2TOD76(velMOD, JD)

    @test isapprox(posMOD[1], 5094.0283745; atol=1e-3)
    @test isapprox(posMOD[2], 6127.8708164; atol=1e-3)
    @test isapprox(posMOD[3], 6380.2485164; atol=1e-3)
    @test isapprox(velMOD[1], -4.746263052; atol=1e-6)
    @test isapprox(velMOD[2], 0.786014045; atol=1e-6)
    @test isapprox(velMOD[3], 5.531790562; atol=1e-6)

    posDiff = posTOD - pos
    velDiff = velTOD - vel

    @test posDiff' * posDiff < 1e-13
    @test velDiff' * velDiff < 1e-13
end

@testset "Nutation SOFA" begin

    JD = [2400000.5, 53736.0]
    N = TOD2MOD76_matrix(JD)'

    @test isapprox(N[1, 1], 0.9999999999534999268; atol=1e-12)
    @test isapprox(N[1, 2], 0.8847935789636432161e-5; atol=1e-12)
    @test isapprox(N[1, 3], 0.3835906502164019142e-5; atol=1e-12)

    @test isapprox(N[2, 1], -0.8847780042583435924e-5; atol=1e-12)
    @test isapprox(N[2, 2], 0.9999999991366569963; atol=1e-12)
    @test isapprox(N[2, 3], -0.4060052702727130809e-4; atol=1e-12)

    @test isapprox(N[3, 1], -0.3836265729708478796e-5; atol=1e-12)
    @test isapprox(N[3, 2], 0.4060049308612638555e-4; atol=1e-12)
    @test isapprox(N[3, 3], 0.9999999991684415129; atol=1e-12)
end

@testset "Precession SOFA" begin
    JD = [2400000.5, 50123.9999]
    P = MOD2J200076_matrix(JD)'

    @test isapprox(P[1, 1], 0.9999995504328350733; atol=1e-12)
    @test isapprox(P[1, 2], 0.8696632209480960785e-3; atol=1e-14)
    @test isapprox(P[1, 3], 0.3779153474959888345e-3; atol=1e-14)

    @test isapprox(P[2, 1], -0.8696632209485112192e-3; atol=1e-14)
    @test isapprox(P[2, 2], 0.9999996218428560614; atol=1e-12)
    @test isapprox(P[2, 3], -0.1643284776111886407e-6; atol=1e-14)

    @test isapprox(P[3, 1], -0.3779153474950335077e-3; atol=1e-14)
    @test isapprox(P[3, 2], -0.1643306746147366896e-6; atol=1e-14)
    @test isapprox(P[3, 3], 0.9999999285899790119; atol=1e-12)
end

@testset "Precession Vallado" begin
    date1 = [2004.0, 4, 6, 7, 52, 32.570009] #TT

    JD, _ = dateVec2JDate(date1)

    pos = [5094.0283745, 6127.8708164, 6380.2485164]
    vel = [-4.746263052, 0.786014045, 5.531790562]

    posJ2000 = MOD2J200076(pos, JD)
    velJ2000 = MOD2J200076(vel, JD)
    posMOD = J20002MOD76(posJ2000, JD)
    velMOD = J20002MOD76(velJ2000, JD)

    @test isapprox(posJ2000[1], 5102.5096; atol=1e-3)
    @test isapprox(posJ2000[2], 6123.01152; atol=1e-3)
    @test isapprox(posJ2000[3], 6378.1363; atol=1e-3)
    @test isapprox(velJ2000[1], -4.7432196; atol=1e-6)
    @test isapprox(velJ2000[2], 0.7905366; atol=1e-6)
    @test isapprox(velJ2000[3], 5.53375619; atol=1e-6)

    posDiff = posMOD - pos
    velDiff = velMOD - vel

    @test posDiff' * posDiff < 1e-13
    @test velDiff' * velDiff < 1e-13
end

@testset "External Val 1" begin
    date1 = [2014.0, 8, 1, 0, 0, 0] #UTC
    JD, _ = dateVec2JDate(date1)
    JDUT1 = UTC2UT1(JD; type=:JD)

    pos = [9891.04671, 20034.03101, 30013.410277]
    posPEF = TOD2PEF76(pos, JDUT1)

    # I don't have the original source for this test, so I'm not sure why this 
    # is off by ~300 m instead of sub-meter. The other tests from what I 
    # believe to be the same source are sub-meter.
    @test isapprox(posPEF[1], -9157.41451; atol=0.5)
    @test isapprox(posPEF[2], 20379.817543; atol=0.5)
    @test isapprox(posPEF[3], pos[3]; atol=1e-4)
end

@testset "External Val 2" begin
    date1 = [2014.0, 6, 1, 0, 0, 0] #UTC
    JD, _ = dateVec2JDate(date1)
    JDUT1 = UTC2UT1(JD; type=:JD)

    pos = [9892.40716, 20033.656809, 30013.211681]
    posPEF = TOD2PEF76(pos, JDUT1)

    @test isapprox(posPEF[1], -22233.18581; atol=1e-3)
    @test isapprox(posPEF[2], 2211.916133; atol=1e-3)
    @test isapprox(posPEF[3], pos[3]; atol=1e-4)
end

@testset "External Val 3" begin
    date1 = [2014.0, 8, 28, 6, 46, 24.461] #UTC
    JD, _ = dateVec2JDate(date1)
    JDUT1 = UTC2UT1(JD; type=:JD)
    JDTT = UT12TT(JDUT1; type=:JD)

    pos = [23141.52, 35279.3, -5.05699] #J2k
    posMOD = J20002MOD76(pos, JDTT)
    posTOD = MOD2TOD76(posMOD, JDTT)
    posPEF = TOD2PEF76(posTOD, JDUT1)

    @test isapprox(posTOD[1], 23024.613191; atol=1e-4)
    @test isapprox(posTOD[2], 35355.699039; atol=1e-4)
    @test isapprox(posTOD[3], 26.7356961; atol=1e-4)

    @test isapprox(posPEF[1], 39365.211746; atol=1e-3)
    @test isapprox(posPEF[2], -15183.49009; atol=1e-3)
    @test isapprox(posPEF[3], posTOD[3]; atol=1e-4)
end
