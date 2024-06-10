# This file uses routines and computations derived by from software provided by
# SOFA under license (see the LICENSE); and does not itself constitute software
# provided by and/or endorsed by SOFA.

"""
    JD, MJD = dateVec2JDate(dateVec)

Convert a date vector into Julian Date and Modified Julian Date.

A date vector is a length 6 vector of floats with
[year, month, day, hour, minute, seconds]

The return values are Julian Date returned in two pieces, in the usual SOFA
manner, which is designed to preserve time resolution. The full Julian Date is
available as a single number by adding the two components of the vector.

Specifically, Julian Date contains the full number of days in JD[1] and the day
fraction in JD[2].

Derived from SOFA's `cal2jd`
"""
function dateVec2JDate(dateVec::Vector{Float64})
    if dateVec[1] < -4799
        error("Year field out of bounds")
    elseif dateVec[2] < 1 || dateVec[2] > 12
        error("Month field out of bounds")
    end
    leapYear =
        dateVec[2] == 2 &&
        dateVec[1] % 4 == 0 &&
        (dateVec[1] % 100 != 0 || dateVec[1] % 400 == 0)
    mtab = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    if dateVec[3] < 1 || dateVec[3] > (mtab[Int(dateVec[2])] + leapYear)
        error("Day field out of bounds")
    end

    my = trunc((dateVec[2] - 14) / 12)
    ypmy = dateVec[1] + my
    mjd =
        trunc(1461 * (ypmy + 4800) / 4) + trunc(367 * (dateVec[2] - 2 - 12 * my) / 12) -
        trunc(3 * ((ypmy + 4900) / 100) / 4) + dateVec[3] - 2432076
    jd = mjd + JM0
    frac = dateVec[4] / 24.0 + dateVec[5] / 1440.0 + dateVec[6] / 86400.0
    return ([jd, frac], [mjd, frac])
end
#TODO: Add a test for this.

"""
    dateVec = JDate2dateVec(JD)

Convert a Julian Date into a date vector

The input values are Julian Date returned in two pieces, in the usual SOFA
manner, which is designed to preserve time resolution. The full Julian Date is
available as a single number by adding the two components of the vector.

The use of the full Julian Date or Modified Julian date is specified by the
`type` option, with `:JD` (default) specifying the full Julian Date, and `:MJD`
specifying Modified Julian Date.

A date vector is a length 6 vector of floats with
[year, month, day, hour, minute, seconds]

Combines SOFA's `jd2cal` and Vallado (v5) Algorithm 22, Vallado for yr, month,
and day (SOFA was giving me day issues), and SOFA for hr, min, sec to preserve
seconds precision.
"""
function JDate2dateVec(JD::Vector{Float64}; type::Symbol=:JD)

    mtab = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    # Convert to full JD if MJD
    useJD = copy(JD)
    if type == :MJD
        useJD[1] += JM0
    end
    j = sum(useJD)

    # Check for date range, algorithm only valid for (1900,2100)
    if j < 2415385.5 || j > 2488068.5
        error("Julian Date out of range, only accepts 1901-2099")
    end

    T1900 = (j - 2415019.5) / 365.25
    yr = 1900 + trunc(T1900)

    leapYrs = trunc((yr - 1900 - 1) * 0.25)
    days = (j - 2415019.5) - ((yr - 1900) * 365.0 + leapYrs)

    if days < 1
        yr -= 1
        leapYrs = trunc((yr - 1900 - 1) * 0.25)
        days = (j - 2415019.5) - ((yr - 1900) * 365.0 + leapYrs)
    end

    if yr % 4 == 0
        mtab[2] = 29
    end

    dayOfYear = trunc(days)
    month = 0
    while sum(mtab[1:(month+1)]) < dayOfYear
        month += 1
    end
    day = dayOfYear - sum(mtab[1:month])
    month += 1

    #Separate date from fraction (-.5 < f < .5)
    d = round(useJD[1])
    f1 = useJD[1] - d
    jd = trunc(d)
    d = round(useJD[2])
    f2 = useJD[2] - d
    jd += trunc(d)

    # Kahan summation for f1+f2+.5 (Klein 2006)
    s = 0.5
    cs = 0.0
    v = [f1, f2]
    for vi in v
        t = s + vi
        cs += abs(s) >= abs(vi) ? (s - t) + vi : (vi - t) + s
        s = t
        if s >= 1.0
            jd += 1
            s -= 1.0
        end
    end
    f = s + cs
    cs = f - s

    # Deal with negative f
    if f < 0
        f = s + 1.0
        cs += (1 - f) + s
        s = f
        f = s + cs
        cs = f - s
        jd -= 1
    end

    # Deal with f that is 1.0 or more when rounded to double
    if f - 1 >= eps(Float64) / 4.0
        t = s - 1.0
        cs += (s - t) - 1.0
        s = t
        f = s + cs
        if -eps(Float64) / 2.0 < f
            jd += 1
            f = max(f, 0.0)
        end
    end
    τ = (f) * 24
    hr = trunc(τ)
    min = trunc((τ - hr) * 60)
    sec = (τ - hr - min / 60) * 3600
    return [yr, month, day, hr, min, sec]
end
# This is the SOFA version, which currently gives incorrect day values.
#function JDate2dateVec(JD::Vector{Float64}; type::Symbol=:JD)
#    # Convert to full JD if MJD
#    j = copy(JD)
#    if type == :MJD
#        j[1] += JM0
#    end
#
#    #Separate date from fraction (-.5 < f < .5)
#    dj = sum(j)
#    d = round(j[1])
#    f1 = j[1] - d
#    jd = trunc( d)
#    d = round(j[2])
#    f2 = j[2] - d
#    jd += trunc( d)
#
#    # Kahan summation for f1+f2+.5 (Klein 2006)
#    s = 0.5;
#    cs = 0.0;
#    v = [f1, f2]
#    for vi in v
#        t = s + vi
#        cs += abs(s) >= abs(vi) ? (s-t) + vi : (vi-t) + s
#        s = t
#        if s >= 1.0
#            jd += 1
#            s -= 1.0
#        end
#    end
#    f = s + cs
#    cs = f - s
#
#    # Deal with negative f
#    if f < 0
#        f = s + 1.
#        cs += (1-f) + s
#        s = f
#        f = s + cs
#        cs = f - s
#        jd -= 1
#    end
#
#    # Deal with f that is 1.0 or more when rounded to double
#    if f-1 >= eps(Float64)/4.
#        t = s - 1.0
#        cs += (s-t) - 1.
#        s = t
#        f = s + cs
#        if -eps(Float64)/2. < f
#            jd += 1
#            f = max(f, 0.0)
#        end
#    end
#
#    # Convert to datevec
#    dateVec = zeros(6)
#    l = jd + 68569
#    n = trunc((4*l)/146097)
#    l -= trunc((146097*n+3)/4)
#    i = trunc((4000 * (l + 1)) / 1461001)
#    l -= trunc((1461*i)/4)-31
#    k = trunc((80*l)/2447)
#    dateVec[3] = trunc((l - (2447*k))/80)
#    l = trunc(k/11)
#    dateVec[2] = k + 2 - 12 * l
#    dateVec[1] = 100 * (n-49) + i + l
#
#    dateVec[4] = trunc(f * 24)
#    f -= dateVec[4]/24
#    dateVec[5] = trunc(f * 1440)
#    f -= dateVec[5]/1440
#    dateVec[6] = f * 86400
#
#    return dateVec
#end

"""
    newDateVec = fixDateVec(dateVec)

Adjust a date vector to account for out of bounds conditions.

For example, [2024., 5, 32, 0, 0, 0] will be adjusted to [2024., 6, 1, 0, 0, 0]

Currently accounts for leap years but not leap seconds. May still have issues
if required to cross multiple leap years.
"""
function fixDateVec(dateVec::Vector{Float64})
    temp = copy(dateVec)
    fixDateVec!(temp)
    return temp
end

"""
    fixDateVec!(dateVec)

Adjust a date vector to account for out of bounds conditions.

For example, [2024., 5, 32, 0, 0, 0] will be adjusted to [2024., 6, 1, 0, 0, 0]

Currently accounts for leap years but not leap seconds. May still have issues
if required to cross multiple leap years.
"""
function fixDateVec!(dateVec::Vector{Float64})

    mtab = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    if dateVec[6] >= 60 || dateVec[6] < 0
        temp = floor(dateVec[6] / 60.0)
        dateVec[5] += temp
        dateVec[6] -= temp * 60.0
    end
    if dateVec[5] >= 60 || dateVec[5] < 0
        temp = floor(dateVec[5] / 60.0)
        dateVec[4] += temp
        dateVec[5] -= temp * 60.0
    end
    if dateVec[4] >= 24 || dateVec[4] < 0
        temp = floor(dateVec[4] / 24.0)
        dateVec[3] += temp
        dateVec[4] -= temp * 24.0
    end

    if dateVec[2] > 12 || dateVec[2] < 1
        temp = floor(dateVec[2] / 12)
        dateVec[1] += temp
        dateVec[2] -= temp * 12.0
    end

    dayOfYear = sum(mtab[1:(Int(dateVec[2])-1)]) + dateVec[3]
    yr = dateVec[1]
    leapYear(x) = x % 4 == 0 && (x % 100 != 0 || x % 400 == 0)

    while dayOfYear < 1
        yr -= 1
        if leapYear(yr)
            dayOfYear += 366.0
        else
            dayOfYear += 365.0
        end
    end
    days = leapYear(yr) ? 366.0 : 365.0
    while dayOfYear > days
        yr += 1
        if leapYear(yr)
            dayOfYear -= days
        else
            dayOfYear -= days
        end
        days = leapYear(yr) ? 366.0 : 365.0
    end

    mtab[2] += leapYear(yr)

    dateVec[1] = yr
    i = 1
    while dayOfYear > mtab[i]
        dayOfYear -= mtab[i]
        i += 1
    end
    dateVec[2] = i
    dateVec[3] = dayOfYear
end
#TODO: test out additional use cases to make sure this is consistent.

"""
    JD_TAI = UTC2TAI(JD_UTC)

Convert a UTC Julian Date into TAI.

The values are Julian Date returned in two pieces, in the usual SOFA
manner, which is designed to preserve time resolution. The full Julian Date is
available as a single number by adding the two components of the vector.

Specifically, Julian Date contains the full number of days in JD[1] and the day
fraction in JD[2].

The `type` optional input specifies if the input `JD` is the full Julian Date
(`:JD`) or the Modified (`:MJD`)

Derived from SOFA's `utctai`
"""
function UTC2TAI(JD::Vector{Float64}; type::Symbol=:JD)
    # Get ΔAT at 0h
    dv = JDate2dateVec(JD; type=type)
    dat0 = DATdateVec(vcat(dv[1:3], zeros(3)))

    # Get ΔAT at 12h
    dat12 = DATdateVec(vcat(dv[1:3], [12.0, 0, 0]))

    # Get ΔAT at 0h tomorrow to detect jumps
    dat24 = DATdateVec(fixDateVec(vcat(dv[1:2], dv[3] + 1, zeros(3))))

    # Separate ΔAT change inter per-day and any jump
    dlod = 2.0 * (dat12 - dat0)
    dleap = dat24 - (dat0 + dlod)

    # Remove any scaling applied to spread leap into preceding day
    fd = dv[6] / 86400 + dv[5] / 1440 + dv[4] / 24
    fd *= (86400 + dleap) / 86400

    # Scale from pre-1972 UTC seconds to SI seconds
    fd *= (86400 + dlod) / 86400

    # Today's calendar date to 2-part JD
    j, m = dateVec2JDate(vcat(dv[1:3], zeros(3)))
    z = type == :JD ? j : m

    # Assemble the TAI result
    tai = zeros(2)
    tai[1] = JD[1]
    a2 = z[1] - JD[1]
    a2 += z[2]
    a2 += fd + dat0 / 86400
    tai[2] = a2

    return tai
end


"""
    JD_TAI = TT2TAI(JD_TT)

Convert a TT Julian Date into TAI.

The values are Julian Date returned in two pieces, in the usual SOFA
manner, which is designed to preserve time resolution. The full Julian Date is
available as a single number by adding the two components of the vector.

Specifically, Julian Date contains the full number of days in JD[1] and the day
fraction in JD[2].

Derived from SOFA's `tttai`
"""
function TT2TAI(JD::Vector{Float64})
    dtat = TTMTAI / 86400.0

    tai = copy(JD)
    tai[2] -= dtat
    return tai
end

"""
    JD_UTC = UT12UTC(JD_UT1)

Convert a UT1 Julian Date into UTC.

The values are Julian Date returned in two pieces, in the usual SOFA
manner, which is designed to preserve time resolution. The full Julian Date is
available as a single number by adding the two components of the vector.

Specifically, Julian Date contains the full number of days in JD[1] and the day
fraction in JD[2].

The use of the full Julian Date or Modified Julian date is specified by the
`type` option, with `:JD` specifying the full Julian Date, and `:MJD` (default)
specifying Modified Julian Date.

Derived from SOFA's `ut1utc`
"""
function UT12UTC(JD::Vector{Float64}; type::Symbol=:MJD)
    u1 = JD[1]
    u2 = JD[2]

    duts = dut1(JD, type=type)
    dats1 = 0.0
    d1 = u1

    for i = -1:3
        d2 = u2 + i
        dateVec = JDate2dateVec(JD, type=type)
        dats2 = DATdateVec(vcat(dateVec[1:3], zeros(3)))
        if i == -1
            dats1 = dats2
        end
        ddats = dats2 - dats1
        if abs(ddats) >= 0.5
            # Leap second nearby: Ensure UT1-UTC is before value
            if ddats * duts >= 0.0
                duts -= ddats
            end

            # UT1 for the start of the UTC day that ends in a leap
            j, m = dateVec2JDate(vcat(dateVec[1:3], zeros(3)))
            if type == :MJD
                d1 = m[1]
                d2 = m[2]
            else
                d1 = j[1]
                d2 = j[2]
            end
            us1 = d1
            us2 = d2 - 1.0 + duts / 86400

            # Is the UT1 after this point?
            du = u1 - us1
            du += u2 - us2
            if du > 0
                # Yes: fraction of the current UTC day that has elapsed
                fd = du * 86400 / (86400 + ddats)
                # Ramp UT1-UTC to bring about SOFA's JD(UTC) convention
                duts += ddats * (fd <= 1.0 ? fd : 1.0)
            end
            break
        end
        dats1 = dats2
    end
    u2 -= duts / 86400.0

    return [u1, u2]
end

"""
    JD_UT1 = UTC2UT1(JD_UTC)

Convert a UTC Julian Date into UT1.

The values are Julian Date returned in two pieces, in the usual SOFA
manner, which is designed to preserve time resolution. The full Julian Date is
available as a single number by adding the two components of the vector.

Specifically, Julian Date contains the full number of days in JD[1] and the day
fraction in JD[2].

The use of the full Julian Date or Modified Julian date is specified by the
`type` option, with `:JD` specifying the full Julian Date, and `:MJD` (default)
specifying Modified Julian Date.

Derived from SOFA's `utcut1`
"""
function UTC2UT1(JD::Vector{Float64}; type::Symbol=:MJD)
    JDTAI = UTC2TAI(JD; type=type)
    return TAI2UT1(JDTAI)
end

"""
    JD_UT1 = TAI2UT1(JD_TAI)

Convert a TAI Julian Date into UT1.

The values are Julian Date returned in two pieces, in the usual SOFA
manner, which is designed to preserve time resolution. The full Julian Date is
available as a single number by adding the two components of the vector.

Specifically, Julian Date contains the full number of days in JD[1] and the day
fraction in JD[2].

The use of the full Julian Date or Modified Julian date is specified by the
`type` option, with `:JD` specifying the full Julian Date, and `:MJD` (default)
specifying Modified Julian Date.

Derived from SOFA's `taiut1`
"""
function TAI2UT1(JD::Vector{Float64}; type::Symbol=:MJD)
    if type == :MJD
        ΔAT = DAT(JD)
    else
        ΔAT = DAT([JD[1] + JM0, JD[2]])
    end
    Δut1 = dut1(JD; type=type)
    dta = Δut1 - ΔAT
    return [JD[1], JD[2] + dta / 86400.0]
end

"""
    JD_TAI = UT12TAI(JD_UT1)

Convert a UT1 Julian Date into TAI.

The values are Julian Date returned in two pieces, in the usual SOFA
manner, which is designed to preserve time resolution. The full Julian Date is
available as a single number by adding the two components of the vector.

Specifically, Julian Date contains the full number of days in JD[1] and the day
fraction in JD[2].

The use of the full Julian Date or Modified Julian date is specified by the
`type` option, with `:JD` specifying the full Julian Date, and `:MJD` (default)
specifying Modified Julian Date.

Derived from SOFA's `ut1tai`
"""
function UT12TAI(JD::Vector{Float64}; type::Symbol=:MJD)
    JD_UTC = UT12UTC(JD, type=type)
    return UTC2TAI(JD_UTC, type=type)
end

"""
    JD_TT = UT12TT(JD_UT1)

Convert a UT1 Julian Date into TT.

The values are Julian Date returned in two pieces, in the usual SOFA
manner, which is designed to preserve time resolution. The full Julian Date is
available as a single number by adding the two components of the vector.

Specifically, Julian Date contains the full number of days in JD[1] and the day
fraction in JD[2].

The use of the full Julian Date or Modified Julian date is specified by the
`type` option, with `:JD` specifying the full Julian Date, and `:MJD` (default)
specifying Modified Julian Date.

Derived from SOFA's `ut1tt`
"""
function UT12TT(JD::Vector{Float64}; type::Symbol=:MJD)
    JD_TAI = UT12TAI(JD, type=type)
    JD_TAI[2] += TTMTAI / 86400.0
    return JD_TAI
end

"""
    JD_UT1 = TT2UT1(JD_TT)

Convert a TT Julian Date into UT1.

The values are Julian Date returned in two pieces, in the usual SOFA
manner, which is designed to preserve time resolution. The full Julian Date is
available as a single number by adding the two components of the vector.

Specifically, Julian Date contains the full number of days in JD[1] and the day
fraction in JD[2].

The use of the full Julian Date or Modified Julian date is specified by the
`type` option, with `:JD` specifying the full Julian Date, and `:MJD` (default)
specifying Modified Julian Date.

Derived from SOFA's `ttut1`
"""
function TT2UT1(JD::Vector{Float64}; type::Symbol=:MJD)
    JD_TAI = TT2TAI(JD)
    return TAI2UT1(JD_TAI; type=type)
end

"""
    JD_TDB = TT2TDB(JD_TT)

Convert a TT Julian Date into TDB.

The values are Julian Date returned in two pieces, in the usual SOFA
manner, which is designed to preserve time resolution. The full Julian Date is
available as a single number by adding the two components of the vector.

Specifically, Julian Date contains the full number of days in JD[1] and the day
fraction in JD[2].

The use of the full Julian Date or Modified Julian date is specified by the
`type` option, with `:JD` (default) specifying the full Julian Date, and `:MJD`
specifying Modified Julian Date.

Derived from SOFA's `tttdb` and 2012 Astronomical Almanac (via Vallado)
"""
function TT2TDB(JD::Vector{Float64}; type::Symbol=:JD)
    JD_TDB = copy(JD)

    Δλ = 246.11 + 0.90251792 * (JD_TDB[1] - 2451545.0 + JD_TDB[2])
    Δλ *= pi / 180
    if type == :MJD
        T = JulianCentury([JD_TDB[1] + JM0, JD_TDB[2]])
    else
        T = JulianCentury(JD_TDB)
    end
    Mearth = 357.5277233 + 35999.05034 * T
    Mearth *= pi / 180

    delta = 0.001657 * sin(Mearth) + 0.000022 * sin(Δλ)
    delta /= 86400.0

    JD_TDB[2] += delta
    return JD_TDB
end

# Requires the full JD, not modified
function JulianCentury(JD::Vector{Float64})
    return (JD[1] - 2451545.0 + JD[2]) / 36525
end
