function wrapTo2pi(input::Float64)
    p2 = 2 * pi
    while input < 0
        input += p2
    end

    return input % p2
end

function wrapToPi(input::Float64)
    p2 = 2 * pi
    w = input % (p2)
    if abs(w) > 0.5
        w -= sign(input) * p2
    end
    return w
end

function R1(theta::Float64)
    R = zeros(3, 3)
    R[1, 1] = 1.0
    ct = cos(theta)
    st = sin(theta)
    R[2, 2] = ct
    R[3, 3] = ct
    R[2, 3] = st
    R[3, 2] = -st

    return R
end

function R2(theta::Float64)
    R = zeros(3, 3)
    R[2, 2] = 1.0
    ct = cos(theta)
    st = sin(theta)
    R[1, 1] = ct
    R[3, 3] = ct
    R[3, 1] = st
    R[1, 3] = -st

    return R
end

function R3(theta::Float64)
    R = zeros(3, 3)
    R[3, 3] = 1.0
    ct = cos(theta)
    st = sin(theta)
    R[2, 2] = ct
    R[1, 1] = ct
    R[1, 2] = st
    R[2, 1] = -st

    return R
end

function cross(a::Vector{Float64}, b::Vector{Float64})
    if !(length(a) == length(b) == 3)
        throw(DimensionMismatch("cross product is only defined for vectors of length 3"))
    end
    a1, a2, a3 = a
    b1, b2, b3 = b

    return [a2 * b3 - a3 * b2, a3 * b1 - a1 * b3, a1 * b2 - a2 * b1]
end
