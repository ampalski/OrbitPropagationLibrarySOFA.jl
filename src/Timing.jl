# This file uses routines and computations derived from software provided by
# SOFA under license (see the LICENSE); and does not itself constitute software
# provided by and/or endorsed by SOFA.

"""
    JD, MJD = datevec2jdate(dateVec; system::Symbol=:UTC)

Convert a date vector into Julian Date and Modified Julian Date.

A date vector is a length 6 vector of floats with
[year, month, day, hour, minute, seconds]

The `system` keyword argument defines which time system the date vector is in.
This can be selected from :UTC (default), :UT1, :TAI, :TT, and :TDB

The return values are Julian Date returned in two pieces, in the usual SOFA
manner, which is designed to preserve time resolution. The full Julian Date is
available as a single number by adding the two components of the vector. This 
can be found in `(M)JD.epoch`

Derived from SOFA's `cal2jd` and `dtf2d`
"""
function datevec2jdate(dateVec::Vector{Float64}; system::Symbol=:UTC)
    if dateVec[1] < -4799
        error("Year field out of bounds")
    elseif dateVec[2] < 1 || dateVec[2] > 12
        error("Month field out of bounds")
    elseif !(system in validSystems)
        error("Invalid time system selection. Valid options are :UTC, :UT1, :TAI, :TT, and :TDB")
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

    # Handle UTC values
    DAYSEC = 86400.0
    seclim = 60.0
    if system == :UTC
        m = MJDate(SA[mjd, 0.0], :UTC)
        dat0 = dat(m)
        m = MJDate(SA[mjd, 0.5], :UTC)
        dat12 = dat(m)
        m = MJDate(SA[mjd+1, 0.0], :UTC)
        dat24 = dat(m)
        dleap = dat24 - (2.0 * dat12 - dat0)
        DAYSEC += dleap
        if dateVec[4] == 23 && dateVec[5] == 59
            seclim += dleap
        end
    end
    frac = (60.0 * (60 * dateVec[4] + dateVec[5]) + dateVec[6]) / DAYSEC

    # frac = dateVec[4] / 24.0 + dateVec[5] / 1440.0 + dateVec[6] / 86400.0
    return (JDate(SA[jd, frac], system), MJDate(SA[mjd, frac], system))
end

"""
    dateVec = jdate2datevec(JD)

Convert a Julian Date into a date vector

The input values are Julian Date returned in two pieces, in the usual SOFA
manner, which is designed to preserve time resolution. The full Julian Date is
available as a single number by adding the two components of the vector.

A date vector is a length 6 vector of floats with
[year, month, day, hour, minute, seconds]

Combines SOFA's `jd2cal` and Vallado (v5) Algorithm 22, Vallado for yr, month,
and day (SOFA was giving me day issues), and SOFA for hr, min, sec to preserve
seconds precision. Also utilized `d2dtf` from SOFA to handle UTC values
"""
function jdate2datevec(JD::JulianDate)

    mtab = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    # Convert to full JD if MJD

    useJD = JD isa JDate ? copy(JD.epoch) : SA[JD.epoch[1]+JM0, JD.epoch[2]]
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

    leap = false
    if JD.system == :UTC
        temp = [yr, month, day, 0.0, 0.0, 0]
        _, m = datevec2jdate(temp, system=:UTC)
        dat0 = dat(m)
        m = MJDate(SA[m.epoch[1], 0.5], m.system)
        dat12 = dat(m)
        m = MJDate(SA[m.epoch[1]+1, 0.0], m.system)
        dat24 = dat(m)
        dleap = dat24 - (2.0 * dat12 - dat0)
        leap = abs(dleap) > 0.5
        if leap
            f += f * dleap / 86400.0
        end
    end

    # Provisional time of day
    τ = (f) * 24
    hr = trunc(τ)
    min = trunc((τ - hr) * 60)
    sec = (τ - hr - min / 60) * 3600
    if hr > 23

        if leap
            if sec >= 1 #past the leap second, go to the next day
                day += 1
                hr = 0.0
                min = 0.0
                sec = 0.0 + (sec - trunc(sec))
                return fixdatevec([yr, month, day, hr, min, sec])
            else
                hr = 23.0
                min = 59.0
                sec = 60.0 + (sec - trunc(sec))
            end
        else
            return fixdatevec([yr, month, day, hr, min, sec])
        end
    end
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
    newDateVec = fixdatevec(dateVec)

Adjust a date vector to account for out of bounds conditions.

For example, [2024., 5, 32, 0, 0, 0] will be adjusted to [2024., 6, 1, 0, 0, 0]

Currently accounts for leap years but not leap seconds. May still have issues
if required to cross multiple leap years.
"""
function fixdatevec(dateVec::Vector{Float64})
    temp = copy(dateVec)
    fixdatevec!(temp)
    return temp
end

"""
    fixdatevec!(dateVec)

Adjust a date vector to account for out of bounds conditions.

For example, [2024., 5, 32, 0, 0, 0] will be adjusted to [2024., 6, 1, 0, 0, 0]

Currently accounts for leap years but not leap seconds. May still have issues
if required to cross multiple leap years.
"""
function fixdatevec!(dateVec::Vector{Float64})

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

    yr = dateVec[1]
    leapYear(x) = x % 4 == 0 && (x % 100 != 0 || x % 400 == 0)
    mtab[2] += leapYear(yr)
    dayOfYear = sum(mtab[1:(Int(dateVec[2])-1)]) + dateVec[3]

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

    dateVec[1] = yr
    i = 1
    while dayOfYear > mtab[i]
        dayOfYear -= mtab[i]
        i += 1
    end
    dateVec[2] = i
    dateVec[3] = dayOfYear
end

function convert_jd(JD::JulianDate, newSystem::Symbol)
    if !(newSystem in validSystems)
        error("Invalid time system selection. Valid options are :UTC, :UT1, :TAI, :TT, and :TDB")
    end
    if JD.system == newSystem
        return JD
    elseif JD.system == :TDB
        error("Conversion from TDB to other systems currently not supported.")
    end

    # utc < > ut1
    # utc > tai
    # ut1 > tt
    # ut1 > tai
    # tai > ut1
    # tt > tai
    # tt > ut1
    # tt > tdb

    # Convert first to UT1, then out to newSystem
    if JD.system == :UT1
        midJD = JD
    elseif JD.system == :UTC
        #UCT2UT1
        midJD = _utc2ut1(JD)
    elseif JD.system == :TAI
        #TAI2UT1
        midJD = _tai2ut1(JD)
    elseif JD.system == :TT
        #TT2UT1
        midJD = _tt2ut1(JD)
    end

    # Convert intermediate UT1 point to newSystem
    if newSystem == :UT1
        return midJD
    elseif newSystem == :UTC
        #UT12UTC
        return _ut12utc(midJD)
    elseif newSystem == :TAI
        #UT12TAI
        return _ut12tai(midJD)
    elseif newSystem == :TT
        # UT12TT
        return _ut12tt(midJD)
    elseif newSystem == :TDB
        # UT12TT -> TT2TDB
        midJD2 = _ut12tt(midJD)
        return _tt2tdb(midJD2)
    else
        error("Undefined output")
    end
end

"""
    JD_TAI = _utc2tai(JD_UTC)

Convert a UTC Julian Date into TAI.

The values are Julian Date returned in two pieces, in the usual SOFA
manner, which is designed to preserve time resolution. The full Julian Date is
available as a single number by adding the two components of the vector.

The `type` optional input specifies if the input `JD` is the full Julian Date
(`:JD`) or the Modified (`:MJD`)

Derived from SOFA's `utctai`
"""
function _utc2tai(JD::JulianDate)
    # Get ΔAT at 0h
    dv = jdate2datevec(JD)
    dat0 = dat_datevec(vcat(dv[1:3], zeros(3)))

    # Get ΔAT at 12h
    dat12 = dat_datevec(vcat(dv[1:3], [12.0, 0, 0]))

    # Get ΔAT at 0h tomorrow to detect jumps
    dat24 = dat_datevec(fixdatevec(vcat(dv[1:2], dv[3] + 1, zeros(3))))

    # Separate ΔAT change inter per-day and any jump
    dlod = 2.0 * (dat12 - dat0)
    dleap = dat24 - (dat0 + dlod)

    # Remove any scaling applied to spread leap into preceding day
    fd = dv[6] / 86400 + dv[5] / 1440 + dv[4] / 24
    fd *= (86400 + dleap) / 86400

    # Scale from pre-1972 UTC seconds to SI seconds
    fd *= (86400 + dlod) / 86400

    # Today's calendar date to 2-part JD
    j, m = datevec2jdate(vcat(dv[1:3], zeros(3)))
    z = JD isa JDate ? j.epoch : m.epoch

    # Assemble the TAI result
    a2 = z[1] - JD.epoch[1]
    a2 += z[2]
    a2 += fd + dat0 / 86400
    if JD isa JDate
        tai = JDate(SA[JD.epoch[1], a2], :TAI)
    else
        tai = MJDate(SA[JD.epoch[1], a2], :TAI)
    end
    return tai
end


"""
    JD_TAI = _tt2tai(JD_TT)

Convert a TT Julian Date into TAI.

The values are Julian Date returned in two pieces, in the usual SOFA
manner, which is designed to preserve time resolution. The full Julian Date is
available as a single number by adding the two components of the vector.

Specifically, Julian Date contains the full number of days in JD[1] and the day
fraction in JD[2].

Derived from SOFA's `tttai`
"""
function _tt2tai(JD::JulianDate)
    dtat = TTMTAI / 86400.0

    if JD isa JDate
        tai = JDate(SA[JD.epoch[1], JD.epoch[2]-dtat], :TAI)
    else
        tai = MJDate(SA[JD.epoch[1], JD.epoch[2]-dtat], :TAI)
    end
    return tai
end

"""
    JD_UTC = _ut12utc(JD_UT1)

Convert a UT1 Julian Date into UTC.

The values are Julian Date returned in two pieces, in the usual SOFA
manner, which is designed to preserve time resolution. The full Julian Date is
available as a single number by adding the two components of the vector.

Derived from SOFA's `ut1utc`
"""
function _ut12utc(JD::JulianDate)
    u1 = JD.epoch[1]
    u2 = JD.epoch[2]

    duts = dut1(JD)
    dats1 = 0.0
    d1 = u1

    for i = -1:3
        d2 = u2 + i
        dateVec = jdate2datevec(JD)
        dats2 = dat_datevec(vcat(dateVec[1:3], zeros(3)))
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
            j, m = datevec2jdate(vcat(dateVec[1:3], zeros(3)))
            if JD isa MJDate
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
    if JD isa JDate
        utc = JDate(SA[u1, u2], :UTC)
    else
        utc = MJDate(SA[u1, u2], :UTC)
    end
    return utc
end

"""
    JD_UT1 = _utc2ut1(JD_UTC)

Convert a UTC Julian Date into UT1.

The values are Julian Date returned in two pieces, in the usual SOFA
manner, which is designed to preserve time resolution. The full Julian Date is
available as a single number by adding the two components of the vector.

Derived from SOFA's `utcut1`
"""
function _utc2ut1(JD::JulianDate)
    JDTAI = _utc2tai(JD)
    return _tai2ut1(JDTAI)
end

"""
    JD_UT1 = _tai2ut1(JD_TAI)

Convert a TAI Julian Date into UT1.

The values are Julian Date returned in two pieces, in the usual SOFA
manner, which is designed to preserve time resolution. The full Julian Date is
available as a single number by adding the two components of the vector.

Derived from SOFA's `taiut1`
"""
function _tai2ut1(JD::JulianDate)
    if JD isa MJDate
        ΔAT = dat(JD)
    else
        ΔAT = dat(MJDate(SA[JD.epoch[1]-JM0, JD.epoch[2]], JD.system))
    end
    Δut1 = dut1(JD)
    dta = Δut1 - ΔAT
    if JD isa JDate
        ut1 = JDate(SA[JD.epoch[1], JD.epoch[2]+dta/86400.0], :UT1)
    else
        ut1 = MJDate(SA[JD.epoch[1], JD.epoch[2]+dta/86400.0], :UT1)
    end
    return ut1
end

"""
    JD_TAI = _ut12tai(JD_UT1)

Convert a UT1 Julian Date into TAI.

The values are Julian Date returned in two pieces, in the usual SOFA
manner, which is designed to preserve time resolution. The full Julian Date is
available as a single number by adding the two components of the vector.

Derived from SOFA's `ut1tai`
"""
function _ut12tai(JD::JulianDate)
    JD_UTC = _ut12utc(JD)
    return _utc2tai(JD_UTC)
end

"""
    JD_TT = _ut12tt(JD_UT1)

Convert a UT1 Julian Date into TT.

The values are Julian Date returned in two pieces, in the usual SOFA
manner, which is designed to preserve time resolution. The full Julian Date is
available as a single number by adding the two components of the vector.

Derived from SOFA's `ut1tt`
"""
function _ut12tt(JD::JulianDate)
    JD_TAI = _ut12tai(JD)
    if JD isa JDate
        tt = JDate(SA[JD_TAI.epoch[1], JD_TAI.epoch[2]+TTMTAI/86400.0], :TT)
    else
        tt = MJDate(SA[JD_TAI.epoch[1], JD_TAI.epoch[2]+TTMTAI/86400.0], :TT)
    end
    return tt
end

"""
    JD_UT1 = _tt2ut1(JD_TT)

Convert a TT Julian Date into UT1.

The values are Julian Date returned in two pieces, in the usual SOFA
manner, which is designed to preserve time resolution. The full Julian Date is
available as a single number by adding the two components of the vector.

Derived from SOFA's `ttut1`
"""
function _tt2ut1(JD::JulianDate)
    JD_TAI = _tt2tai(JD)
    return _tai2ut1(JD_TAI)
end

"""
    JD_TDB = _tt2tdb(JD_TT)

Convert a TT Julian Date into TDB.

The values are Julian Date returned in two pieces, in the usual SOFA
manner, which is designed to preserve time resolution. The full Julian Date is
available as a single number by adding the two components of the vector.

Derived from SOFA's `tttdb` and 2012 Astronomical Almanac (via Vallado)
"""
function _tt2tdb(JD::JulianDate)
    if JD isa JDate
        temp = JD.epoch[1] - 2451545.0 + JD.epoch[2]
        T = juliancentury(JD)
    else
        temp = JD.epoch[1] + JM0 - 2451545.0 + JD.epoch[2]
        T = juliancentury(mjdate_to_jdate(JD))
    end

    Δλ = 246.11 + 0.90251792 * temp
    Δλ *= pi / 180
    Mearth = 357.5277233 + 35999.05034 * T
    Mearth *= pi / 180

    delta = 0.001657 * sin(Mearth) + 0.000022 * sin(Δλ)
    delta /= 86400.0

    if JD isa JDate
        tdb = JDate(SA[JD.epoch[1], JD.epoch[2]+delta], :TDB)
    else
        tdb = MJDate(SA[JD.epoch[1], JD.epoch[2]+delta], :TDB)
    end
    return tdb
end

# Requires the full JD, not modified
function juliancentury(JD::JDate)
    return (JD.epoch[1] - J00 + JD.epoch[2]) / 36525
end

"""
    GMST = gmst(JD_UT1)

Convert a UT1 Julian Date into the GMST at that epoch.

The input values are Julian Date returned in two pieces, in the usual SOFA
manner, which is designed to preserve time resolution. The full Julian Date is
available as a single number by adding the two components of the vector.

The Greenwich Mean Sidereal Time is an angle describing the mean rotation of
the Prime Meridian from the vernal equinox. There are several implementations,
which can be selected by the `model` input, 1982 (`82`, default), 2000 (`00`)
or 2006(`06`).

Derived from SOFA's `gmst82`, `gmst00`, and `gmst06`
"""
function gmst(JD::JulianDate; model::Integer=82)
    useJD = JD isa JDate ? JD : mjdate_to_jdate(JD)
    useJD = useJD.system == :UT1 ? useJD : convert_jd(useJD, :UT1)

    if model == 82
        return _gmst82(useJD)
    end
    JDtt = _ut12tt(useJD)
    if model == 00
        return _gmst00(useJD, JDtt)
    else
        return _gmst06(useJD, JDtt)
    end
end

function _gmst82(JD::JDate)
    A = 24110.54841 - 43200.0
    B = 8640184.812866
    C = 0.093104
    D = -6.2e-6

    t = juliancentury(JD)
    f = 86400 * (JD.epoch[1] % 1.0 + JD.epoch[2] % 1.0)

    return wrapto2pi(S2R * ((A + (B + (C + D * t) * t) * t) + f))
end
function _gmst00(JD::JDate, JDtt::JDate)
    t = juliancentury(JDtt)
    era = _era00(JD)
    return wrapto2pi(era +
                     (0.014506 +
                      (4612.15739966 +
                       (1.39667721 -
                        (0.00009344 +
                         (0.00001882) * t) * t) * t) * t) * AS2R)
end
function _gmst06(JD::JDate, JDtt::JDate)
    t = juliancentury(JDtt)
    era = _era00(JD)
    return wrapto2pi(era +
                     (0.014506 +
                      (4612.156534 +
                       (1.3915817 -
                        (0.00000044 -
                         (0.000029956 -
                          (0.0000000368) * t) * t) * t) * t) * t) * AS2R)
end

function _era00(JD::JDate)
    t = JD.epoch[2] + (JD.epoch[1] - 2451545.0)

    f = JD.epoch[1] % 1.0 + JD.epoch[2] % 1.0

    theta = wrapto2pi(2 * pi * (f + 0.7790572732640 + 0.00273781191135448 * t))

    return theta
end

"""
    eqeq = eqeq94(JD)

Convert a TT Julian Date into the 1994 Equation of the Equinoxes

The input values are Julian Date returned in two pieces, in the usual SOFA
manner, which is designed to preserve time resolution. The full Julian Date is
available as a single number by adding the two components of the vector.

Derived from SOFA's `eqeq94`
"""
function eqeq94(JD::JulianDate)
    useJD = JD isa JDate ? JD : jdate_to_mjdate(JD)
    om = _lunarlan(useJD)
    dψ, _ = nutationterms(useJD)
    eps0 = obl(useJD; model=80)
    ee = dψ * cos(eps0) + AS2R * (0.00264 * sin(om) + 0.000063 * sin(om + om))

    return ee
end


"""
    GAST = gast(JD_UT1)

Convert a UT1 Julian Date into the GAST at that epoch.

The input values are Julian Date returned in two pieces, in the usual SOFA
manner, which is designed to preserve time resolution. The full Julian Date is
available as a single number by adding the two components of the vector.

The Greenwich Apparent Sidereal Time is an angle describing the true rotation 
of the Prime Meridian from the vernal equinox. There are several
implementations, which can be selected by the `model` input, 1982/1994 (`94`, 
default), 2000 (`00`) or 2006(`06`).

Derived from SOFA's `gst94`, `gmst00`, and `gmst06`
"""
function gast(JD::JulianDate; model::Integer=94)
    JDUT1 = JD isa JDate ? JD : mjdate_to_jdate(JD)
    JDTT = _ut12tt(JDUT1)

    if model == 94
        GMST = gmst(JDUT1, model=82)
        eqeq = eqeq94(JDTT)
        return wrapto2pi(GMST + eqeq)
    end
end



