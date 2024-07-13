# This file uses routines and computations derived from software provided by
# SOFA under license (see the LICENSE); and does not itself constitute software
# provided by and/or endorsed by SOFA.

module OrbitPropagationLibrarySOFA

using DataFrames
using JLD2, FileIO
using StaticArrays

export JulianDate, JDate, MJDate
export jdate_to_mjdate, mjdate_to_jdate
export dateVec2JDate, JDate2dateVec
export fixDateVec, fixDateVec!
export DAT, DATdateVec
export convert_jd
# export _utc2tai, _tt2tai, _tai2ut1
# export _ut12utc, _utc2ut1
# export _ut12tai, _ut12tt
# export _tt2ut1, _tt2tdb
export JulianCentury
export gmst, gast
export ITRF2PEF76_matrix, ITRF2PEF76, PEF2ITRF76
export PEF2TOD76_matrix, PEF2TOD76, PEF2TOD76_vel, TOD2PEF76, TOD2PEF76_vel
export TOD2MOD76_matrix, TOD2MOD76, MOD2TOD76
export MOD2J200076_matrix, MOD2J200076, J20002MOD76
export TEME2TOD_matrix, TEME2TOD, TOD2TEME
# export dut1, EOP, JM0 # only use for debugging
# export obl, NutationTerms, eqeq94 #use only for debugging
# Include constituent files
include("Utils.jl")
include("TypeDefs.jl")
include("Constants.jl")
include("Timing.jl")
include("Coordinates.jl")

end
