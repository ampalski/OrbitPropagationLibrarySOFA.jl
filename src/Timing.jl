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

    dayOfYear = sum(mtab[1:(dateVec[2] - 1)]) + dateVec[3]
    yr = dateVec[1]
    leapYear(x) = x % 4 == 0 && (x % 100 != 0 || x % 400 == 0)

    while dayOfYear < 1
        yr -= 1
        if leapYear(yr)
            dayOfYear += 366.
        else
            dayOfYear += 365.
        end
    end
    days = leapYear(yr) ? 366. : 365.
    while dayOfYear > days 
        yr += 1
        if leapYear(yr)
            dayOfYear -= days
        else
            dayOfYear -= days
        end
        days = leapYear(yr) ? 366. : 365.
    end

    mtab[2] += leapYear(yr)
    
    dateVec[1] = yr
    i = 1
    while dayOfYear > mtab[i]
        dayOfYear -= mtab[i]
        i+=1
    end
    dateVec[2] = i
    dateVec[3] = dayOfYear
end
#TODO: test out additional use cases to make sure this is consistent.
