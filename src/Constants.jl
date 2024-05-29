# This file uses routines and computations derived by from software provided by
# SOFA under license (see the LICENSE); and does not itself constitute software
# provided by and/or endorsed by SOFA.

# add to this from sofam.h as needed

const WGS84 = 1
const GRS80 = 2
const WGS72 = 3

const R2AS = 206264.8062470963551564734
const AS2R = 4.848136811095359935899141e-6

const J00 = 2451545.0
const JM0 = 2400000.5
const JM00 = 51544.5

const TTMTAI = 32.184

const AU = 149597870.7 #km
const CMPS = 299792458.0 #m/s

"""
    ΔAT = dat(dateVec)

Find the TAI-UTC value (ΔAT) from a date vector

A date vector is a length 6 vector of floats with
[year, month, day, hour, minute, seconds]

    TODO: Finish this description. based on iauDat.c
The return values are Julian Date returned in two pieces, in the usual SOFA
manner, which is designed to preserve time resolution. The full Julian Date is
available as a single number by adding the two components of the vector.

Specifically, Julian Date contains the full number of days in JD[1] and the day
fraction in JD[2].

Derived from SOFA's `cal2jd`
"""
function dat(dateVec::Vector{Float64})

end

