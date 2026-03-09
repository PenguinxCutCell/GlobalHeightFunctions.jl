spatial_shape_from_grid(grid::CartesianGrids.CartesianGrid{N}) where {N} = grid.n

function reshape_V(V::AbstractVector{T}, grid::CartesianGrids.CartesianGrid{N}) where {T,N}
    length(V) == prod(grid.n) || throw(DimensionMismatch("vector length $(length(V)) must match prod(grid.n)=$(prod(grid.n))"))
    return reshape(copy(V), grid.n...)
end

function column_profile(V::AbstractArray{T,N}; axis::Union{Int,Symbol}=1) where {T,N}
    ia = axis_to_index(axis, Val(N))
    if N == 1
        out = Array{T,0}(undef)
        out[] = sum(V)
        return out
    end
    return dropdims(sum(V, dims=ia), dims=ia)
end

function _transverse_measure(grid::CartesianGrids.CartesianGrid{N,T}, ia::Int) where {N,T}
    Δ = CartesianGrids.meshsize(grid)
    m = one(T)
    @inbounds for d in 1:N
        d == ia && continue
        m *= convert(T, Δ[d])
    end
    return m
end

function positions_from_heights(H, grid::CartesianGrids.CartesianGrid{N,T}; axis::Union{Int,Symbol}=1) where {N,T}
    ia = axis_to_index(axis, Val(N))
    x0 = convert(T, first(CartesianGrids.grid1d(grid, ia)))
    scale = _transverse_measure(grid, ia)
    return x0 .+ H ./ scale
end

function ensure_periodic!(xf::AbstractVector)
    isempty(xf) && return xf
    xf[end] = xf[1]
    return xf
end

function ensure_periodic!(xf::AbstractMatrix)
    ny, nz = size(xf)
    if ny > 1
        xf[end, :] .= xf[1, :]
    end
    if nz > 1
        xf[:, end] .= xf[:, 1]
    end
    return xf
end

ensure_periodic!(xf) = xf
