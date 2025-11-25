# This file uses routines and computations derived from software provided by
# SOFA under license (see the LICENSE); and does not itself constitute software
# provided by and/or endorsed by SOFA.

const validFrames = [:ITRF, :PEF, :TOD, :TEME, :MOD, :J2000]
function _num_frame(frame::Symbol)
    if frame == :ITRF
        return 1
    elseif frame == :PEF
        return 2
    elseif frame == :TOD || frame == :TEME
        return 3
    elseif frame == :MOD
        return 4
    else
        return 5
    end
end

"""
    pos_new = convert_pos(pos, oldFrame, newFrame, epoch)

Convert an orbital position given in `oldFrame` to `newFrame` at a given `epoch`

Allowable frames are: `:ITRF`, `:PEF`, `:TOD`, `:TEME`, `:MOD`, `:J2000`

The input value for `epoch` is the Julian Date returned in two pieces, in the 
usual SOFA manner, which is designed to preserve time resolution. The full 
Julian Date is available as a single number by adding the two components of the
vector.

# Inputs
* `pos::AbstractVector`: 3-element position vector for the given velocity, in `oldFrame`
* `oldFrame::Symbol`: The frame for the input position and velocity
* `newFrame::Symbol`: The output frame to convert the velocity to
* `epoch::JulianDate`: The Julian Date for the specified position and velocity

"""
function convert_pos(
    pos::AbstractVector,
    oldFrame::Symbol,
    newFrame::Symbol,
    epoch::JulianDate
)
    if !(oldFrame in validFrames)
        error("Invalid source coordinate frame provided.")
    end
    oldFrameNum = _num_frame(oldFrame)
    if !(newFrame in validFrames)
        error("Invalid destination coordinate frame provided.")
    end
    newFrameNum = _num_frame(newFrame)
    if newFrame == oldFrame
        return pos
    end
    if pos isa SVector{3,Float64}
        svFlag = true
        usePos = pos[1:3]
    else
        svFlag = false
        if length(pos) != 3
            error("Position vectors must be of length 3")
        end
        usePos = copy(pos)
    end

    if oldFrame == :TEME
        epochuse = convert_jd(epoch, :TT)
        usePos = teme2tod(usePos, epochuse)
    end

    currentFrame = oldFrameNum
    dir = newFrameNum > oldFrameNum
    ktr = 0
    while currentFrame != newFrameNum
        if currentFrame == 1
            epochuse = convert_jd(epoch, :UTC)
            usePos = itrf2pef76(usePos, epochuse)
            currentFrame = 2
        elseif currentFrame == 2 && dir
            epochuse = convert_jd(epoch, :UT1)
            usePos = pef2tod76(usePos, epochuse)
            currentFrame = 3
        elseif currentFrame == 2
            epochuse = convert_jd(epoch, :UTC)
            usePos = pef2itrf76(usePos, epochuse)
            currentFrame = 1
        elseif currentFrame == 3 && dir
            epochuse = convert_jd(epoch, :TT)
            usePos = tod2mod76(usePos, epochuse)
            currentFrame = 4
        elseif currentFrame == 3
            epochuse = convert_jd(epoch, :UT1)
            usePos = tod2pef76(usePos, epochuse)
            currentFrame = 2
        elseif currentFrame == 4 && dir
            epochuse = convert_jd(epoch, :TT)
            usePos = mod2j200076(usePos, epochuse)
            currentFrame = 5
        elseif currentFrame == 4
            epochuse = convert_jd(epoch, :TT)
            usePos = mod2tod76(usePos, epochuse)
            currentFrame = 3
        else
            epochuse = convert_jd(epoch, :TT)
            usePos = j20002mod76(usePos, epochuse)
            currentFrame = 4
        end
        ktr += 1
        if ktr > 8
            error("Inf looped")
        end
    end

    if newFrame == :TEME
        epochuse = convert_jd(epoch, :TT)
        usePos = tod2teme(usePos, epochuse)
    end

    if svFlag
        return SA[usePos...]
    else
        return usePos
    end
end

"""
    (pos_new, vel_new) = convert_posvel(pos, vel, oldFrame, newFrame, epoch)

Convert an orbital position and velocity given in `oldFrame` to `newFrame` at a given `epoch`

Allowable frames are: `:ITRF`, `:PEF`, `:TOD`, `:TEME`, `:MOD`, `:J2000`

The input value for `epoch` is the Julian Date returned in two pieces, in the 
usual SOFA manner, which is designed to preserve time resolution. The full 
Julian Date is available as a single number by adding the two components of the
vector.

# Inputs
* `pos::AbstractVector`: 3-element position vector for the given velocity, in `oldFrame`
* `vel::AbstractVector`: 3-element velocity vector, in `oldFrame`
* `oldFrame::Symbol`: The frame for the input position and velocity
* `newFrame::Symbol`: The output frame to convert the velocity to
* `epoch::JulianDate`: The Julian Date for the specified position and velocity

"""
function convert_posvel(
    pos::AbstractVector,
    vel::AbstractVector,
    oldFrame::Symbol,
    newFrame::Symbol,
    epoch::JulianDate
)
    if !(oldFrame in validFrames)
        error("Invalid source coordinate frame provided.")
    end
    oldFrameNum = _num_frame(oldFrame)
    if !(newFrame in validFrames)
        error("Invalid destination coordinate frame provided.")
    end
    newFrameNum = _num_frame(newFrame)
    if newFrame == oldFrame
        return (pos, vel)
    end
    if pos isa SVector{3,Float64} && vel isa SVector{3,Float64}
        svFlag = true
        usePos = pos[1:3]
        useVel = vel[1:3]
    else
        svFlag = false
        if length(pos) != 3
            error("Position vectors must be of length 3")
        end
        if length(vel) != 3
            error("Velocity vectors must be of length 3")
        end
        usePos = copy(pos)
        useVel = copy(vel)
    end

    if oldFrame == :TEME
        epochuse = convert_jd(epoch, :TT)
        usePos = teme2tod(usePos, epochuse)
        useVel = teme2tod(useVel, epochuse)
    end

    currentFrame = oldFrameNum
    dir = newFrameNum > oldFrameNum
    ktr = 0
    while currentFrame != newFrameNum
        if currentFrame == 1
            epochuse = convert_jd(epoch, :UTC)
            usePos = itrf2pef76(usePos, epochuse)
            useVel = itrf2pef76(useVel, epochuse)
            currentFrame = 2
        elseif currentFrame == 2 && dir
            epochuse = convert_jd(epoch, :UT1)
            useVel = pef2tod76_vel(useVel, usePos, epochuse)
            usePos = pef2tod76(usePos, epochuse)
            currentFrame = 3
        elseif currentFrame == 2
            epochuse = convert_jd(epoch, :UTC)
            usePos = pef2itrf76(usePos, epochuse)
            useVel = pef2itrf76(useVel, epochuse)
            currentFrame = 1
        elseif currentFrame == 3 && dir
            epochuse = convert_jd(epoch, :TT)
            usePos = tod2mod76(usePos, epochuse)
            useVel = tod2mod76(useVel, epochuse)
            currentFrame = 4
        elseif currentFrame == 3
            epochuse = convert_jd(epoch, :UT1)
            usePos = tod2pef76(usePos, epochuse)
            useVel = tod2pef76_vel(useVel, usePos, epochuse)
            currentFrame = 2
        elseif currentFrame == 4 && dir
            epochuse = convert_jd(epoch, :TT)
            usePos = mod2j200076(usePos, epochuse)
            useVel = mod2j200076(useVel, epochuse)
            currentFrame = 5
        elseif currentFrame == 4
            epochuse = convert_jd(epoch, :TT)
            usePos = mod2tod76(usePos, epochuse)
            useVel = mod2tod76(useVel, epochuse)
            currentFrame = 3
        else
            epochuse = convert_jd(epoch, :TT)
            usePos = j20002mod76(usePos, epochuse)
            useVel = j20002mod76(useVel, epochuse)
            currentFrame = 4
        end
        ktr += 1
        if ktr > 8
            error("Inf looped")
        end
    end

    if newFrame == :TEME
        epochuse = convert_jd(epoch, :TT)
        usePos = tod2teme(usePos, epochuse)
        useVel = tod2teme(useVel, epochuse)
    end

    if svFlag
        return (SA[usePos...], SA[useVel...])
    else
        return (usePos, useVel)
    end
end

"""
    state_new = convert_state(state, oldFrame, newFrame, epoch)

Convert an orbital state given in `oldFrame` to `newFrame` at a given `epoch`

Allowable frames are: `:ITRF`, `:PEF`, `:TOD`, `:TEME`, `:MOD`, `:J2000`

The input value for `epoch` is the Julian Date returned in two pieces, in the 
usual SOFA manner, which is designed to preserve time resolution. The full 
Julian Date is available as a single number by adding the two components of the
vector.

# Inputs
* `state::AbstractVector`: 6-element position/velocity vector, in `oldFrame`
* `oldFrame::Symbol`: The frame for the input position and velocity
* `newFrame::Symbol`: The output frame to convert the velocity to
* `epoch::JulianDate`: The Julian Date for the specified position and velocity

"""
function convert_state(
    state::AbstractVector,
    oldFrame::Symbol,
    newFrame::Symbol,
    epoch::JulianDate,
)
    if state isa SVector{6,Float64}
        svFlag = true
    else
        svFlag = false
        if length(state) != 6
            error("State vectors must be of length 6")
        end
    end

    usePos = state[1:3]
    useVel = state[4:6]

    usePos, useVel = convert_posvel(usePos, useVel, oldFrame, newFrame, epoch)
    state = vcat(usePos, useVel)
    if svFlag
        return SA[state...]
    else
        return state
    end
end

"""
    vel_new = convert_vel(pos, vel, oldFrame, newFrame, epoch)

Convert a velocity given in `oldFrame` to `newFrame` at a given `epoch`

Allowable frames are: `:ITRF`, `:PEF`, `:TOD`, `:TEME`, `:MOD`, `:J2000`

The input value for `epoch` is the Julian Date returned in two pieces, in the 
usual SOFA manner, which is designed to preserve time resolution. The full 
Julian Date is available as a single number by adding the two components of the
vector.

# Inputs
* `pos::AbstractVector`: 3-element position vector for the given velocity, in `oldFrame`
* `vel::AbstractVector`: 3-element velocity vector, in `oldFrame`
* `oldFrame::Symbol`: The frame for the input position and velocity
* `newFrame::Symbol`: The output frame to convert the velocity to
* `epoch::JulianDate`: The Julian Date for the specified position and velocity

"""
function convert_vel(
    pos::AbstractVector,
    vel::AbstractVector,
    oldFrame::Symbol,
    newFrame::Symbol,
    epoch::JulianDate
)
    _, newVel = convert_posvel(pos, vel, oldFrame, newFrame, epoch)
    return newVel
end
