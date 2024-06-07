# This file uses routines and computations derived by from software provided by
# SOFA under license (see the LICENSE); and does not itself constitute software
# provided by and/or endorsed by SOFA.

# add to this from sofam.h as needed
EOP = load("data/EOP.jld2")["EOP"]

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
    dat = DAT(MJD)

Return the ΔAT value for the given MJD.

The input MJD values are Modified Julian Date returned in two pieces, in the
usual SOFA manner, which is designed to preserve time resolution. The full 
Julian Date is available as a single number by adding the two components of the
vector. For example, [2024., 5, 32, 0, 0, 0] will be adjusted to 
[2024., 6, 1, 0, 0, 0]

Requires an update with every new leap second.

Derived from SOFA's `iauDat`
"""
function DAT(MJD::Vector{Float64})
    drift = [
        37300.0 0.0012960;
        37300.0 0.0012960;
        37300.0 0.0012960;
        37665.0 0.0011232;
        37665.0 0.0011232;
        38761.0 0.0012960;
        38761.0 0.0012960;
        38761.0 0.0012960;
        38761.0 0.0012960;
        38761.0 0.0012960;
        38761.0 0.0012960;
        38761.0 0.0012960;
        39126.0 0.0025920;
        39126.0 0.002592
    ]
    changes = [
        36934.0 1.417818;
        37300.0 1.422818;
        37512.0 1.372818;
        37665.0 1.845858;
        38334.0 1.945858;
        38395.0 3.240130;
        38486.0 3.340130;
        38639.0 3.440130;
        38761.0 3.540130;
        38820.0 3.640130;
        38942.0 3.740130;
        39004.0 3.840130;
        39126.0 4.313170;
        39887.0 4.213170;
        41317.0 10.0;
        41499.0 11.0;
        41683.0 12.0;
        42048.0 13.0;
        42413.0 14.0;
        42778.0 15.0;
        43144.0 16.0;
        43509.0 17.0;
        43874.0 18.0;
        44239.0 19.0;
        44786.0 20.0;
        45151.0 21.0;
        45516.0 22.0;
        46247.0 23.0;
        47161.0 24.0;
        47892.0 25.0;
        48257.0 26.0;
        48804.0 27.0;
        49169.0 28.0;
        49534.0 29.0;
        50083.0 30.0;
        50630.0 31.0;
        51179.0 32.0;
        53736.0 33.0;
        54832.0 34.0;
        56109.0 35.0;
        57204.0 36.0;
        57754.0 37.0
    ]
    MJDbase = floor(sum(MJD))

    # if pre-UTC year, error out
    if MJDbase < changes[1, 1]
        error("Pre-UTC year")
    end

    index = findlast(MJDbase .> @view changes[:, 1])

    if isnothing(index)
        error("Underflow, no data available")
    end

    dat = changes[index, 2]

    #if pre-1972, adjust for drift
    if MJDbase < 41317
        dat += (sum(MJD) - drift[index, 1]) * drift[index, 2]
    end

    return dat
end


"""
    dat = DATdateVec(dateVec)

Return the ΔAT value for the given date vector.

A date vector is a length 6 vector of floats with
[year, month, day, hour, minute, seconds]

Requires an update with every new leap second.

Derived from SOFA's `iauDat`
"""
function DATdateVec(dateVec::Vector{Float64})
    _, MJD = dateVec2JDate(dateVec)
    return DAT(MJD)
end

"""
    dut1 = dut1(MJD)

Return the UTC-UT1 value for the given MJD date.

The input values are Julian Date returned in two pieces, in the usual SOFA
manner, which is designed to preserve time resolution. The full Julian Date is
available as a single number by adding the two components of the vector.

The use of the full Julian Date or Modified Julian date is specified by the
`type` option, with `:JD` specifying the full Julian Date, and `:MJD` (default)
specifying Modified Julian Date.
"""
function dut1(MJD::Vector{Float64}; type::Symbol=:MJD)
    if type != :MJD
        MJD[1] -= JM0
    end

    firstDate = EOP[1, :MJD] - 1
    date = Int(floor(MJD[1] + MJD[2])) - firstDate

    return EOP[date, :dUT1]
end
