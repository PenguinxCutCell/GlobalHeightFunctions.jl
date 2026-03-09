@inline function _clamp_coord(v::T, a::T, b::T) where {T}
    return clamp(v, min(a, b), max(a, b))
end

@inline function _wrap_periodic(v::T, a::T, b::T) where {T}
    L = b - a
    if L <= zero(T)
        return v
    end
    return a + mod(v - a, L)
end

function _interp_pwlinear_1d(nodes::AbstractVector{T}, vals::AbstractVector{T}, s::T) where {T}
    n = length(nodes)
    n == length(vals) || throw(DimensionMismatch("interpolation nodes and values must have same length"))
    n == 1 && return vals[1]

    if s <= nodes[1]
        return vals[1]
    elseif s >= nodes[end]
        return vals[end]
    end

    i = clamp(searchsortedlast(nodes, s), 1, n - 1)
    s0 = nodes[i]
    s1 = nodes[i + 1]
    θ = (s - s0) / (s1 - s0)
    return (one(T) - θ) * vals[i] + θ * vals[i + 1]
end

function _interp_bilinear(
    ynodes::AbstractVector{T},
    znodes::AbstractVector{T},
    vals::AbstractMatrix{T},
    y::T,
    z::T,
) where {T}
    ny = length(ynodes)
    nz = length(znodes)
    ny == size(vals, 1) || throw(DimensionMismatch("values first dimension must match y nodes"))
    nz == size(vals, 2) || throw(DimensionMismatch("values second dimension must match z nodes"))
    ny == 1 && nz == 1 && return vals[1, 1]

    y_eval = _clamp_coord(y, ynodes[1], ynodes[end])
    z_eval = _clamp_coord(z, znodes[1], znodes[end])

    iy = ny == 1 ? 1 : clamp(searchsortedlast(ynodes, y_eval), 1, ny - 1)
    iz = nz == 1 ? 1 : clamp(searchsortedlast(znodes, z_eval), 1, nz - 1)

    y1 = ynodes[iy]
    y2 = ynodes[min(iy + 1, ny)]
    z1 = znodes[iz]
    z2 = znodes[min(iz + 1, nz)]

    q11 = vals[iy, iz]
    q12 = vals[iy, min(iz + 1, nz)]
    q21 = vals[min(iy + 1, ny), iz]
    q22 = vals[min(iy + 1, ny), min(iz + 1, nz)]

    ty = y2 == y1 ? zero(T) : (y_eval - y1) / (y2 - y1)
    tz = z2 == z1 ? zero(T) : (z_eval - z1) / (z2 - z1)

    return (one(T) - ty) * (one(T) - tz) * q11 +
           ty * (one(T) - tz) * q21 +
           (one(T) - ty) * tz * q12 +
           ty * tz * q22
end
