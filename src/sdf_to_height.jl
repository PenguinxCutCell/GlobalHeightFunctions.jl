function _zero_crossing_location(xnodes::AbstractVector{T}, line::AbstractVector{T}) where {T}
    n = length(xnodes)
    n == length(line) || throw(DimensionMismatch("line and node coordinates must have the same length"))

    @inbounds if iszero(line[1])
        return xnodes[1]
    end

    @inbounds for i in 1:(n - 1)
        a = line[i]
        b = line[i + 1]
        if iszero(a)
            return xnodes[i]
        elseif iszero(b)
            return xnodes[i + 1]
        elseif a * b < zero(T)
            θ = a / (a - b)
            return xnodes[i] + θ * (xnodes[i + 1] - xnodes[i])
        end
    end

    all_pos = true
    all_neg = true
    best_i = 1
    best_abs = abs(line[1])
    @inbounds for i in 1:n
        v = line[i]
        all_pos &= v >= zero(T)
        all_neg &= v <= zero(T)
        av = abs(v)
        if av < best_abs
            best_abs = av
            best_i = i
        end
    end

    if all_pos
        return xnodes[1]
    elseif all_neg
        return xnodes[end]
    end
    return xnodes[best_i]
end

function xf_from_sdf(
    phi::AbstractArray{T,N},
    grid::CartesianGrids.CartesianGrid{N,T};
    axis::Union{Int,Symbol}=1,
    method::Symbol=:zero_crossing,
) where {T,N}
    method === :zero_crossing || throw(ArgumentError("unsupported method `$method`; expected :zero_crossing"))
    size(phi) == grid.n || throw(DimensionMismatch("phi size must match grid.n"))

    ia = axis_to_index(axis, Val(N))
    xnodes = collect(T, CartesianGrids.grid1d(grid, ia))

    if N == 1
        out = Array{T,0}(undef)
        out[] = _zero_crossing_location(xnodes, vec(phi))
        return out
    end

    tshape = ntuple(k -> grid.n[k < ia ? k : k + 1], N - 1)
    xf = Array{T,N-1}(undef, tshape...)
    line = Vector{T}(undef, grid.n[ia])

    @inbounds for itrans in CartesianIndices(tshape)
        for ii in 1:grid.n[ia]
            line[ii] = phi[_line_index(ia, itrans, ii, Val(N))]
        end
        xf[itrans] = _zero_crossing_location(xnodes, line)
    end

    return xf
end
