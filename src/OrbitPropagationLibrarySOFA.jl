# This file uses routines and computations derived by from software provided by
# SOFA under license (see the LICENSE); and does not itself constitute software
# provided by and/or endorsed by SOFA.

module OrbitPropagationLibrarySOFA

using DataFrames, JLD2, FileIO

export dateVec2JDate, JDate2dateVec
export fixDateVec, fixDateVec!
export DAT, DATdateVec
export UTC2TAI, TT2TAI, TAI2UT1
export UT12UTC, UTC2UT1
export UT12TAI, UT12TT
export TT2UT1, TT2TDB
export JulianCentury
export GMST, GAST
# export dut1, EOP, JM0 # only use for debugging
# export OBL, NutationTerms #use only for debugging
# Include constituent files
include("Utils.jl")
include("Constants.jl")
include("Timing.jl")
include("Coordinates.jl")

end
