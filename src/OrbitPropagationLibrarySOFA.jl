# This file uses routines and computations derived from software provided by
# SOFA under license (see the LICENSE); and does not itself constitute software
# provided by and/or endorsed by SOFA.

module OrbitPropagationLibrarySOFA

using FileIO, JLD2
using StaticArrays

export JDate, JulianDate, MJDate
export jdate_to_mjdate, mjdate_to_jdate
export datevec2jdate, jdate2datevec
export fixdatevec
export dat, dat_datevec
export convert_jd
export juliancentury
export gast, gmst
export itrf2pef76, itrf2pef76_matrix, pef2itrf76
export pef2tod76, pef2tod76_matrix, pef2tod76_vel, tod2pef76, tod2pef76_vel
export mod2tod76, tod2mod76, tod2mod76_matrix
export j20002mod76, mod2j200076, mod2j200076_matrix
export teme2tod, teme2tod_matrix, tod2teme
export convert_pos, convert_posvel, convert_state, convert_vel

# Include constituent files
include("Utils.jl")
include("TypeDefs.jl")
include("Constants.jl")
include("Timing.jl")
include("Coordinates.jl")
include("CoordinateConversions.jl")

end
