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
