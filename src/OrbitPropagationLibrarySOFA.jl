# This file uses routines and computations derived from software provided by
# SOFA under license (see the LICENSE); and does not itself constitute software
# provided by and/or endorsed by SOFA.

module OrbitPropagationLibrarySOFA

using DataFrames
using JLD2, FileIO
using StaticArrays

export JulianDate, JDate, MJDate
export jdate_to_mjdate, mjdate_to_jdate
export datevec2jdate, jdate2datevec
export fixdatevec, fixdatevec!
export dat, dat_datevec
export convert_jd
# export _utc2tai, _tt2tai, _tai2ut1
# export _ut12utc, _utc2ut1
# export _ut12tai, _ut12tt
# export _tt2ut1, _tt2tdb
export juliancentury
export gmst, gast
export itrf2pef76_matrix, itrf2pef76, pef2itrf76
export pef2tod76_matrix, pef2tod76, pef2tod76_vel, tod2pef76, tod2pef76_vel
export tod2mod76_matrix, tod2mod76, mod2tod76
export mod2j200076_matrix, mod2j200076, j20002mod76
export teme2tod_matrix, teme2tod, tod2teme
# export dut1, EOP, JM0 # only use for debugging
# export obl, nutationterms, eqeq94 #use only for debugging
# Include constituent files
include("Utils.jl")
include("TypeDefs.jl")
include("Constants.jl")
include("Timing.jl")
include("Coordinates.jl")

end
