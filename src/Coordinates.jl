# This file uses routines and computations derived from software provided by
# SOFA under license (see the LICENSE); and does not itself constitute software
# provided by and/or endorsed by SOFA.

"""
    obl = OBL(JD_TT, type::Symbol=:JD, model::Integer=80)

Convert a TT Julian Date into the Mean Obliquity of the Ecliptic

The input values are Julian Date returned in two pieces, in the usual SOFA
manner, which is designed to preserve time resolution. The full Julian Date is
available as a single number by adding the two components of the vector.

Specifically, Julian Date contains the full number of days in JD[1] and the day
fraction in JD[2]. Requires the full Julian Date (type=:JD).

The Mean Obliquity of the Ecliptic has two implementations,which can be 
selected by the `model` input, 1980 (`:80`, default), or 2006 (`:06`)

Derived from SOFA's `obl80` and `obl06`
"""
function OBL(JD::Vector{Float64}; type::Symbol=:JD, model::Integer=80)
    useJD = copy(JD)
    if type != :JD
        useJD[1] += JM0
    end

    if model == 80
        return _obl80(useJD)
    else
        return _obl06(useJD)
    end
end

function _obl80(JD::Vector{Float64})
    t = JulianCentury(JD)

    return AS2R * (84381.448 + (-46.815 + (-0.00059 + 0.001813 * t) * t) * t)
end
function _obl06(JD::Vector{Float64})
    t = JulianCentury(JD)
    eps0 = (84381.406 +
            (-46.836769 +
             (-0.0001831 +
              (0.0020034 +
               (-0.000000576 +
                (-0.0000000434) * t) * t) * t) * t) * t) * AS2R
    return eps0
end

#JD must be :JD, not :MJD, and in TT
function _lunarLAN(JD::Vector{Float64})
    t = JulianCentury(JD)
    om = wrapToPi((450160.28 + (-482890.539 +
                                (7.455 + 0.008 * t) * t) * t) * AS2R +
                  (-5.0 * t) % 1.0 * 2 * pi)
    return om
end


"""
    dψ, dϵ = NutationTerms(JD; type::Symbol=:JD, numTerms::Integer=106)

Convert a TT Julian Date into the nutation terms from the 1980 model.

The input values are Julian Date returned in two pieces, in the usual SOFA
manner, which is designed to preserve time resolution. The full Julian Date is
available as a single number by adding the two components of the vector.

Specifically, Julian Date contains the full number of days in JD[1] and the day
fraction in JD[2]. 

The `numTerms` optional argument governs how many coefficients are used in the
summation series. The default 106 is the full definition, consistent with the
1980 model. TEME and SGP4 use a 4-term approximation.

Derived from SOFA's `nut80`
"""
function NutationTerms(
    JD::Vector{Float64};
    type::Symbol=:JD,
    numTerms::Integer=106,
)

    useJD = copy(JD)
    if type != :JD
        useJD[1] += JM0
    end

    U2R = AS2R / 1e4 # milliarcseconds conversion
    nutTerms = [
        0.0 0.0 0.0 0.0 1.0 -171996.0 -174.2 92025.0 8.9;
        0.0 0.0 2.0 -2.0 2.0 -13187.0 -1.6 5736.0 -3.1;
        0.0 0.0 2.0 0.0 2.0 -2274.0 -0.2 977.0 -0.5;
        0.0 0.0 0.0 0.0 2.0 2062.0 0.2 -895.0 0.5;
        0.0 -1.0 0.0 0.0 0.0 -1426.0 3.4 54.0 -0.1;
        1.0 0.0 0.0 0.0 0.0 712.0 0.1 -7.0 0.0;
        0.0 1.0 2.0 -2.0 2.0 -517.0 1.2 224.0 -0.6;
        0.0 0.0 2.0 0.0 1.0 -386.0 -0.4 200.0 0.0;
        1.0 0.0 2.0 0.0 2.0 -301.0 0.0 129.0 -0.1;
        0.0 -1.0 2.0 -2.0 2.0 217.0 -0.5 -95.0 0.3;
        -1.0 0.0 0.0 2.0 0.0 158.0 0.0 -1.0 0.0;
        0.0 0.0 2.0 -2.0 1.0 129.0 0.1 -70.0 0.0;
        -1.0 0.0 2.0 0.0 2.0 123.0 0.0 -53.0 0.0;
        1.0 0.0 0.0 0.0 1.0 63.0 0.1 -33.0 0.0;
        0.0 0.0 0.0 2.0 0.0 63.0 0.0 -2.0 0.0;
        -1.0 0.0 2.0 2.0 2.0 -59.0 0.0 26.0 0.0;
        -1.0 0.0 0.0 0.0 1.0 -58.0 -0.1 32.0 0.0;
        1.0 0.0 2.0 0.0 1.0 -51.0 0.0 27.0 0.0;
        -2.0 0.0 0.0 2.0 0.0 -48.0 0.0 1.0 0.0;
        -2.0 0.0 2.0 0.0 1.0 46.0 0.0 -24.0 0.0;
        0.0 0.0 2.0 2.0 2.0 -38.0 0.0 16.0 0.0;
        2.0 0.0 2.0 0.0 2.0 -31.0 0.0 13.0 0.0;
        2.0 0.0 0.0 0.0 0.0 29.0 0.0 -1.0 0.0;
        1.0 0.0 2.0 -2.0 2.0 29.0 0.0 -12.0 0.0;
        0.0 0.0 2.0 0.0 0.0 26.0 0.0 -1.0 0.0;
        0.0 0.0 2.0 -2.0 0.0 -22.0 0.0 0.0 0.0;
        -1.0 0.0 2.0 0.0 1.0 21.0 0.0 -10.0 0.0;
        0.0 2.0 0.0 0.0 0.0 17.0 -0.1 0.0 0.0;
        0.0 2.0 2.0 -2.0 2.0 -16.0 0.1 7.0 0.0;
        -1.0 0.0 0.0 2.0 1.0 16.0 0.0 -8.0 0.0;
        0.0 1.0 0.0 0.0 1.0 -15.0 0.0 9.0 0.0;
        1.0 0.0 0.0 -2.0 1.0 -13.0 0.0 7.0 0.0;
        0.0 -1.0 0.0 0.0 1.0 -12.0 0.0 6.0 0.0;
        2.0 0.0 -2.0 0.0 0.0 11.0 0.0 0.0 0.0;
        -1.0 0.0 2.0 2.0 1.0 -10.0 0.0 5.0 0.0;
        1.0 0.0 2.0 2.0 2.0 -8.0 0.0 3.0 0.0;
        0.0 -1.0 2.0 0.0 2.0 -7.0 0.0 3.0 0.0;
        0.0 0.0 2.0 2.0 1.0 -7.0 0.0 3.0 0.0;
        1.0 1.0 0.0 -2.0 0.0 -7.0 0.0 0.0 0.0;
        0.0 1.0 2.0 0.0 2.0 7.0 0.0 -3.0 0.0;
        -2.0 0.0 0.0 2.0 1.0 -6.0 0.0 3.0 0.0;
        0.0 0.0 0.0 2.0 1.0 -6.0 0.0 3.0 0.0;
        2.0 0.0 2.0 -2.0 2.0 6.0 0.0 -3.0 0.0;
        1.0 0.0 0.0 2.0 0.0 6.0 0.0 0.0 0.0;
        1.0 0.0 2.0 -2.0 1.0 6.0 0.0 -3.0 0.0;
        0.0 0.0 0.0 -2.0 1.0 -5.0 0.0 3.0 0.0;
        0.0 -1.0 2.0 -2.0 1.0 -5.0 0.0 3.0 0.0;
        2.0 0.0 2.0 0.0 1.0 -5.0 0.0 3.0 0.0;
        1.0 -1.0 0.0 0.0 0.0 5.0 0.0 0.0 0.0;
        1.0 0.0 0.0 -1.0 0.0 -4.0 0.0 0.0 0.0;
        0.0 0.0 0.0 1.0 0.0 -4.0 0.0 0.0 0.0;
        0.0 1.0 0.0 -2.0 0.0 -4.0 0.0 0.0 0.0;
        1.0 0.0 -2.0 0.0 0.0 4.0 0.0 0.0 0.0;
        2.0 0.0 0.0 -2.0 1.0 4.0 0.0 -2.0 0.0;
        0.0 1.0 2.0 -2.0 1.0 4.0 0.0 -2.0 0.0;
        1.0 1.0 0.0 0.0 0.0 -3.0 0.0 0.0 0.0;
        1.0 -1.0 0.0 -1.0 0.0 -3.0 0.0 0.0 0.0;
        -1.0 -1.0 2.0 2.0 2.0 -3.0 0.0 1.0 0.0;
        0.0 -1.0 2.0 2.0 2.0 -3.0 0.0 1.0 0.0;
        1.0 -1.0 2.0 0.0 2.0 -3.0 0.0 1.0 0.0;
        3.0 0.0 2.0 0.0 2.0 -3.0 0.0 1.0 0.0;
        -2.0 0.0 2.0 0.0 2.0 -3.0 0.0 1.0 0.0;
        1.0 0.0 2.0 0.0 0.0 3.0 0.0 0.0 0.0;
        -1.0 0.0 2.0 4.0 2.0 -2.0 0.0 1.0 0.0;
        1.0 0.0 0.0 0.0 2.0 -2.0 0.0 1.0 0.0;
        -1.0 0.0 2.0 -2.0 1.0 -2.0 0.0 1.0 0.0;
        0.0 -2.0 2.0 -2.0 1.0 -2.0 0.0 1.0 0.0;
        -2.0 0.0 0.0 0.0 1.0 -2.0 0.0 1.0 0.0;
        2.0 0.0 0.0 0.0 1.0 2.0 0.0 -1.0 0.0;
        3.0 0.0 0.0 0.0 0.0 2.0 0.0 0.0 0.0;
        1.0 1.0 2.0 0.0 2.0 2.0 0.0 -1.0 0.0;
        0.0 0.0 2.0 1.0 2.0 2.0 0.0 -1.0 0.0;
        1.0 0.0 0.0 2.0 1.0 -1.0 0.0 0.0 0.0;
        1.0 0.0 2.0 2.0 1.0 -1.0 0.0 1.0 0.0;
        1.0 1.0 0.0 -2.0 1.0 -1.0 0.0 0.0 0.0;
        0.0 1.0 0.0 2.0 0.0 -1.0 0.0 0.0 0.0;
        0.0 1.0 2.0 -2.0 0.0 -1.0 0.0 0.0 0.0;
        0.0 1.0 -2.0 2.0 0.0 -1.0 0.0 0.0 0.0;
        1.0 0.0 -2.0 2.0 0.0 -1.0 0.0 0.0 0.0;
        1.0 0.0 -2.0 -2.0 0.0 -1.0 0.0 0.0 0.0;
        1.0 0.0 2.0 -2.0 0.0 -1.0 0.0 0.0 0.0;
        1.0 0.0 0.0 -4.0 0.0 -1.0 0.0 0.0 0.0;
        2.0 0.0 0.0 -4.0 0.0 -1.0 0.0 0.0 0.0;
        0.0 0.0 2.0 4.0 2.0 -1.0 0.0 0.0 0.0;
        0.0 0.0 2.0 -1.0 2.0 -1.0 0.0 0.0 0.0;
        -2.0 0.0 2.0 4.0 2.0 -1.0 0.0 1.0 0.0;
        2.0 0.0 2.0 2.0 2.0 -1.0 0.0 0.0 0.0;
        0.0 -1.0 2.0 0.0 1.0 -1.0 0.0 0.0 0.0;
        0.0 0.0 -2.0 0.0 1.0 -1.0 0.0 0.0 0.0;
        0.0 0.0 4.0 -2.0 2.0 1.0 0.0 0.0 0.0;
        0.0 1.0 0.0 0.0 2.0 1.0 0.0 0.0 0.0;
        1.0 1.0 2.0 -2.0 2.0 1.0 0.0 -1.0 0.0;
        3.0 0.0 2.0 -2.0 2.0 1.0 0.0 0.0 0.0;
        -2.0 0.0 2.0 2.0 2.0 1.0 0.0 -1.0 0.0;
        -1.0 0.0 0.0 0.0 2.0 1.0 0.0 -1.0 0.0;
        0.0 0.0 -2.0 2.0 1.0 1.0 0.0 0.0 0.0;
        0.0 1.0 2.0 0.0 1.0 1.0 0.0 0.0 0.0;
        -1.0 0.0 4.0 0.0 2.0 1.0 0.0 0.0 0.0;
        2.0 1.0 0.0 -2.0 0.0 1.0 0.0 0.0 0.0;
        2.0 0.0 0.0 2.0 0.0 1.0 0.0 0.0 0.0;
        2.0 0.0 2.0 -2.0 1.0 1.0 0.0 -1.0 0.0;
        2.0 0.0 -2.0 0.0 1.0 1.0 0.0 0.0 0.0;
        1.0 -1.0 0.0 -2.0 0.0 1.0 0.0 0.0 0.0;
        -1.0 0.0 0.0 1.0 1.0 1.0 0.0 0.0 0.0;
        -1.0 -1.0 0.0 2.0 1.0 1.0 0.0 0.0 0.0;
        0.0 1.0 0.0 1.0 0.0 1.0 0.0 0.0 0.0
    ]
    t = JulianCentury(useJD)

    # Fundamental Arguments

    # Mean longitude of Moon minus mean longitude of Moon's perigee
    el = wrapToPi((485866.733 +
                   (715922.633 +
                    (31.310 + 0.064 * t) * t) * t) * AS2R +
                  ((1325 * t) % 1.0) * 2 * pi)

    # Mean longitude of Sun minus mean longitude of Sun's perigee
    elp = wrapToPi((1287099.804 +
                    (1292581.224 +
                     (-0.577 - 0.012 * t) * t) * t) * AS2R +
                   ((99.0 * t) % 1.0) * 2 * pi)

    # Mean longitude of Moon minus mean longitude of Moon's node
    f = wrapToPi((335778.877 +
                  (295263.137 +
                   (-13.257 + 0.011 * t) * t) * t) * AS2R +
                 ((1342.0 * t) % 1.0) * 2 * pi)

    # Mean elongation of Moon from Sun
    d = wrapToPi((1072261.307 +
                  (1105601.328 +
                   (-6.891 + 0.019 * t) * t) * t) * AS2R +
                 ((1236.0 * t) % 1.0) * 2 * pi)

    # Longitude of the mean ascending node of the lunar orbit on the ecliptic,
    # measured from the mean equinox of date
    om = _lunarLAN(useJD)

    # Nutation series
    dp = 0.0
    de = 0.0

    for ii in numTerms:-1:1
        arg = nutTerms[ii, 1] * el +
              nutTerms[ii, 2] * elp +
              nutTerms[ii, 3] * f +
              nutTerms[ii, 4] * d +
              nutTerms[ii, 5] * om

        s = nutTerms[ii, 6] + nutTerms[ii, 7] * t
        c = nutTerms[ii, 8] + nutTerms[ii, 9] * t
        s != 0 && (dp += s * sin(arg))
        c != 0 && (de += c * cos(arg))
    end

    dpsi = dp * U2R
    deps = de * U2R

    return (dpsi, deps)
end

"""
    ζ, z, θ = PrecessionTerms(JD; type::Symbol=:JD, JD0=[J00, 0.0])

Convert a TT Julian Date into the precession terms from the 1976 model.

The input values are Julian Date returned in two pieces, in the usual SOFA
manner, which is designed to preserve time resolution. The full Julian Date is
available as a single number by adding the two components of the vector.

Specifically, Julian Date contains the full number of days in JD[1] and the day
fraction in JD[2]. 

The `JD0` optional input value specifies the start epoch of the rotation, 
typically set to the J2000 epoch.

Derived from SOFA's `prec76`
"""
function PrecessionTerms(
    JD::Vector{Float64};
    type::Symbol=:JD,
    JD0::Vector{Float64}=[J00, 0.0],
)
    # Interval between the fundamental epoch J2000 and the start date.
    t0 = JulianCentury(JD0)

    # Interval over which precession required
    t = ((JD[1] - JD0[1]) + (JD[2] - JD0[2])) / 36525.0

    # Euler Angles
    tas2r = t * AS2R
    w = 2306.2181 + (1.39656 - 0.000139 * t0) * t0
    ζ = (w + ((0.30188 - 0.000344 * t0) + 0.017998 * t) * t) * tas2r
    z = (w + ((1.09468 + 0.000066 * t0) + 0.018203 * t) * t) * tas2r
    θ = ((2004.3109 + (-0.8533 - 0.000217 * t0) * t0) +
         ((-0.42665 - 0.000217 * t0) - 0.041833 * t) * t) * tas2r

    return (ζ, z, θ)
end
"""
    rotMatrix = ITRF2PEF76_matrix(JD)

Find the ITRF to PEF rotation matrix for a given Julian Date using the IAU-76
model.

The input values are Julian Date is given in two pieces, in the usual SOFA
manner, which is designed to preserve time resolution. The full Julian Date is
available as a single number by adding the two components of the vector.

Used to transform an ITRF vector into PEF as `r_PEF = rotMatrix * r_ITRF`

Derived from Vallado's description of the IAU-76 reduction.
"""
function ITRF2PEF76_matrix(JD::Vector{Float64})
    useJD = [JD[1] - JM0, JD[2]]
    firstDate = EOP[1, :MJD] - 1
    date = Int(floor(useJD[1] + useJD[2])) - firstDate

    xp = EOP[date, :xp] * AS2R
    yp = EOP[date, :yp] * AS2R
    cx = cos(xp)
    sx = sin(xp)
    cy = cos(yp)
    sy = sin(yp)
    W = zeros(3, 3)
    W[1, 1] = cx
    W[1, 3] = -sx

    W[2, 1] = sx * sy
    W[2, 2] = cy
    W[2, 3] = cx * sy

    W[3, 1] = sx * cy
    W[3, 2] = -sy
    W[3, 3] = cx * cy

    # W[1, 1] = 1.0 # approximation version
    # W[2, 2] = 1.0
    # W[3, 3] = 1.0
    #
    # W[1, 3] = -xp
    # W[2, 3] = yp
    # W[3, 1] = xp
    # W[3, 2] = -yp

    return W
end

"""
    r_PEF = ITRF2PEF76(r_ITRF, JD)

Transform an ITRF vector into PEF at the given Julian Date using the IAU-76 
model.

The state vector must be of length 3.

The input Julian Date is given in two pieces, in the usual SOFA
manner, which is designed to preserve time resolution. The full Julian Date is
available as a single number by adding the two components of the vector.

Derived from Vallado's description of the IAU-76 reduction.
"""
function ITRF2PEF76(vec::Vector{Float64}, JD::Vector{Float64})
    W = ITRF2PEF76_matrix(JD)
    return W * vec
end

"""
    r_PEF = PEF2ITRF76(r_ITRF, JD)

Transform a PEF vector into ITRF at the given Julian Date using the IAU-76 
model.

The state vector must be of length 3.

The input Julian Date is given in two pieces, in the usual SOFA
manner, which is designed to preserve time resolution. The full Julian Date is
available as a single number by adding the two components of the vector.

Derived from Vallado's description of the IAU-76 reduction.
"""
function PEF2ITRF76(vec::Vector{Float64}, JD::Vector{Float64})
    W = ITRF2PEF76_matrix(JD)
    return W' * vec
end

"""
    rotMatrix = PEF2TOD76_matrix(JD)

Find the PEF to TOD rotation matrix for a given UT1 Julian Date using the
IAU-76 model.

The input values are Julian Date is given in two pieces, in the usual SOFA
manner, which is designed to preserve time resolution. The full Julian Date is
available as a single number by adding the two components of the vector.

Used to transform a PEF vector into TOD as `r_TOD = rotMatrix * r_PEF`

Derived from Vallado's description of the IAU-76 reduction.
"""
function PEF2TOD76_matrix(JD::Vector{Float64})
    gast = GAST(JD, type=:JD, model=94)
    R = R3(-gast)
    return R
end

"""
    r_TOD = PEF2TOD76(r_PEF, JD)

Transform a PEF vector into TOD at the given UT1 Julian Date using the IAU-76 
model.

The state vector must be of length 3.

The input Julian Date is given in two pieces, in the usual SOFA
manner, which is designed to preserve time resolution. The full Julian Date is
available as a single number by adding the two components of the vector.

Derived from Vallado's description of the IAU-76 reduction.
"""
function PEF2TOD76(vec::Vector{Float64}, JD::Vector{Float64})
    R = PEF2TOD76_matrix(JD)
    return R * vec
end

"""
    r_PEF = TOD2PEF76(r_TOD, JD)

Transform a TOD vector into PEF at the given UT1 Julian Date using the IAU-76 
model.

The state vector must be of length 3.

The input Julian Date is given in two pieces, in the usual SOFA
manner, which is designed to preserve time resolution. The full Julian Date is
available as a single number by adding the two components of the vector.

Derived from Vallado's description of the IAU-76 reduction.
"""
function TOD2PEF76(vec::Vector{Float64}, JD::Vector{Float64})
    R = PEF2TOD76_matrix(JD)
    return R' * vec
end

"""
    v_TOD = PEF2TOD76_vel(v_PEF, r_PEF JD)

Transform a PEF velocity vector into TOD at the given Julian Date using the
IAU-76 model.

The state vector must be of length 3.

The input Julian Date is given in two pieces, in the usual SOFA
manner, which is designed to preserve time resolution. The full Julian Date is
available as a single number by adding the two components of the vector.

Derived from Vallado's description of the IAU-76 reduction.
"""
function PEF2TOD76_vel(
    vec::Vector{Float64},
    pos::Vector{Float64},
    JD::Vector{Float64},
)
    R = PEF2TOD76_matrix(JD)
    useJD = [JD[1] - JM0, JD[2]]
    firstDate = EOP[1, :MJD] - 1
    date = Int(floor(useJD[1] + useJD[2])) - firstDate

    LOD = EOP[date, :LOD]

    omega = [0, 0, 7.292115146706979e-5 * (1 - LOD / 86400)]
    return (R * (vec + cross(omega, pos)))
end

"""
    v_PEF = TOD2PEF76_vel(v_TOD, r_PEF, JD)

Transform a TOD velocity vector into PEF at the given Julian Date using the 
IAU-76 model.

The state vector must be of length 3.

The input Julian Date is given in two pieces, in the usual SOFA
manner, which is designed to preserve time resolution. The full Julian Date is
available as a single number by adding the two components of the vector.

Derived from Vallado's description of the IAU-76 reduction.
"""
function TOD2PEF76_vel(
    vec::Vector{Float64},
    pos::Vector{Float64},
    JD::Vector{Float64},
)
    R = PEF2TOD76_matrix(JD)
    useJD = [JD[1] - JM0, JD[2]]
    firstDate = EOP[1, :MJD] - 1
    date = Int(floor(useJD[1] + useJD[2])) - firstDate

    LOD = EOP[date, :LOD]

    omega = [0, 0, 7.292115146706979e-5 * (1 - LOD / 86400)]

    return (R' * vec - cross(omega, pos))
end

"""
    rotMatrix = TOD2MOD76_matrix(JD)

Find the TOD to MOD rotation matrix for a given TT Julian Date using the IAU-76
model.

The input values are Julian Date is given in two pieces, in the usual SOFA
manner, which is designed to preserve time resolution. The full Julian Date is
available as a single number by adding the two components of the vector.

Used to transform a TOD vector into MOD as `r_MOD = rotMatrix * r_TOD`

Derived from SOFA's `nutm80` and Vallado's description of the IAU-76 reduction.
"""
function TOD2MOD76_matrix(JD::Vector{Float64})
    dψ, dϵ = NutationTerms(JD; type=:JD, numTerms=106)
    obl = OBL(JD; type=:JD, model=80)
    N = R1(-obl) * R3(dψ) * R1((obl + dϵ))
    return N
end

"""
    r_MOD = TOD2MOD76(r_TOD, JD)

Transform a TOD vector into MOD at the given TT Julian Date using the IAU-76 
model.

The state vector must be of length 3.

The input Julian Date is given in two pieces, in the usual SOFA
manner, which is designed to preserve time resolution. The full Julian Date is
available as a single number by adding the two components of the vector.

Derived from SOFA's `nutm80` and Vallado's description of the IAU-76 reduction.
"""
function TOD2MOD76(vec::Vector{Float64}, JD::Vector{Float64})
    N = TOD2MOD76_matrix(JD)
    return N * vec
end

"""
    r_TOD = MOD2TOD76(r_MOD, JD)

Transform a MOD vector into TOD at the given TT Julian Date using the IAU-76 
model.

The state vector must be of length 3.

The input Julian Date is given in two pieces, in the usual SOFA
manner, which is designed to preserve time resolution. The full Julian Date is
available as a single number by adding the two components of the vector.

Derived from SOFA's `nutm80` and Vallado's description of the IAU-76 reduction.
"""
function MOD2TOD76(vec::Vector{Float64}, JD::Vector{Float64})
    N = TOD2MOD76_matrix(JD)
    return N' * vec
end

"""
    rotMatrix = MOD2J200076_matrix(JD)

Find the MOD to J2000 rotation matrix for a given TT Julian Date using the
IAU-76 model.

The input values are Julian Date is given in two pieces, in the usual SOFA
manner, which is designed to preserve time resolution. The full Julian Date is
available as a single number by adding the two components of the vector.

Used to transform a MOD vector into J2000 as `r_J2000 = rotMatrix * r_MOD`

Derived from SOFA's `pmat76` and Vallado's description of the IAU-76 reduction.
"""
function MOD2J200076_matrix(JD::Vector{Float64})
    ζ, z, θ = PrecessionTerms(JD; type=:JD)
    P = R3(ζ) * R2(-θ) * R3(z)
    return P
end

"""
    r_J2000 = MOD2J200076(r_MOD, JD)

Transform a MOD vector into J2000 at the given TT Julian Date using the IAU-76 
model.

The state vector must be of length 3.

The input Julian Date is given in two pieces, in the usual SOFA
manner, which is designed to preserve time resolution. The full Julian Date is
available as a single number by adding the two components of the vector.

Derived from SOFA's `pmat76` and Vallado's description of the IAU-76 reduction.
"""
function MOD2J200076(vec::Vector{Float64}, JD::Vector{Float64})
    P = MOD2J200076_matrix(JD)
    return P * vec
end

"""
    r_MOD = J20002MOD76(r_J2000, JD)

Transform a J2000 vector into MOD at the given TT Julian Date using the IAU-76 
model.

The state vector must be of length 3.

The input Julian Date is given in two pieces, in the usual SOFA
manner, which is designed to preserve time resolution. The full Julian Date is
available as a single number by adding the two components of the vector.

Derived from SOFA's `pmat76` and Vallado's description of the IAU-76 reduction.
"""
function J20002MOD76(vec::Vector{Float64}, JD::Vector{Float64})
    P = MOD2J200076_matrix(JD)
    return P' * vec
end

"""
    rotMatrix = TEME2TOD_matrix(JD)

Find the TEME to TOD rotation matrix for a given TT Julian Date using the
IAU-76 model.

The input values are Julian Date is given in two pieces, in the usual SOFA
manner, which is designed to preserve time resolution. The full Julian Date is
available as a single number by adding the two components of the vector.

Used to transform a TEME vector into TOD as `r_TOD = rotMatrix * r_TEME`

Derived from Vallado's description of the TEME frame. Note that TEME is not 
fully defined publicly, and may carry some errors due to the ambiguous nature 
of its definition.
"""
function TEME2TOD_matrix(JD::Vector{Float64})
    eqeq = eqeq94(JD)
    T = R3(-eqeq)
    return T
end

"""
    r_TOD = TEME2TOD(r_TEME, JD)

Transform a TEME vector into TOD at the given TT Julian Date using the IAU-76 
model.

The state vector must be of length 3.

The input Julian Date is given in two pieces, in the usual SOFA
manner, which is designed to preserve time resolution. The full Julian Date is
available as a single number by adding the two components of the vector.

Derived from Vallado's description of the TEME frame. Note that TEME is not 
fully defined publicly, and may carry some errors due to the ambiguous nature 
of its definition.
"""
function TEME2TOD(vec::Vector{Float64}, JD::Vector{Float64})
    T = TEME2TOD_matrix(JD)
    return T * vec
end

"""
    r_TEME = TOD2TEME(r_TOD, JD)

Transform a TOD vector into TEME at the given TT Julian Date using the IAU-76 
model.

The state vector must be of length 3.

The input Julian Date is given in two pieces, in the usual SOFA
manner, which is designed to preserve time resolution. The full Julian Date is
available as a single number by adding the two components of the vector.

Derived from Vallado's description of the TEME frame. Note that TEME is not 
fully defined publicly, and may carry some errors due to the ambiguous nature 
of its definition.
"""
function TOD2TEME(vec::Vector{Float64}, JD::Vector{Float64})
    T = TEME2TOD_matrix(JD)
    return T' * vec
end
