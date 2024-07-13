abstract type JulianDate end

struct JDate <: JulianDate
    epoch::SVector{2,AbstractFloat}
    system::Symbol
end

struct MJDate <: JulianDate
    epoch::SVector{2,AbstractFloat}
    system::Symbol
end

Base.isless(a::JDate, b::JDate) = (a.epoch[1] + a.epoch[2]) < (b.epoch[1] + b.epoch[2])
Base.isless(a::MJDate, b::MJDate) = (a.epoch[1] + a.epoch[2]) < (b.epoch[1] + b.epoch[2])
function Base.isless(a::JDate, b::MJDate)
    return (a.epoch[1] + a.epoch[2]) < (b.epoch[1] + b.epoch[2] + JM0)
end
function Base.isless(a::MJDate, b::JDate)
    return (a.epoch[1] + a.epoch[2] + JM0) < (b.epoch[1] + b.epoch[2])
end

jdate_to_mjdate(j::JDate) = MJDate(SA[j.epoch[1]-JM0, j.epoch[2]], j.system)
mjdate_to_jdate(j::MJDate) = JDate(SA[j.epoch[1]+JM0, j.epoch[2]], j.system)

const validSystems = [:UT1, :UTC, :TAI, :TT, :TDB]

