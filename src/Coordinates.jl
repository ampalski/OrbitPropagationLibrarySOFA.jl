# This file uses routines and computations derived by from software provided by
# SOFA under license (see the LICENSE); and does not itself constitute software
# provided by and/or endorsed by SOFA.

"""
    obl = OBL(JD_TT, model=:80)

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
function OBL(JD::Vector{Float64}; type::Symbol=:JD, model::Symbol=:80)
    useJD = copy(JD)
    if type != :JD
        useJD[1] += JM0
    end

    if model == :80
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
    dψ, dϵ = NutationTerms(JD)

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
    t = JulianCentury(JD)


end

