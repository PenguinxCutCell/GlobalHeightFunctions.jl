function phi_from_xf!(
    phi::AbstractArray{T,N},
    xf,
    grid::CartesianGrids.CartesianGrid{N,T};
    axis::Union{Int,Symbol}=1,
    interp::Symbol=:linear,
    periodic::Bool=false,
) where {T,N}
    size(phi) == grid.n || throw(DimensionMismatch("phi size must match grid.n"))

    ia = axis_to_index(axis, Val(N))
    tshape = N == 1 ? () : ntuple(k -> grid.n[k < ia ? k : k + 1], N - 1)

    if N == 1
        xf0 = xf isa AbstractArray ? xf[] : convert(T, xf)
        xnodes = collect(T, CartesianGrids.grid1d(grid, ia))
        @inbounds for i in eachindex(phi)
            phi[i] = xnodes[i] - xf0
        end
        return phi
    end

    xf_arr = xf isa AbstractArray ? xf : reshape([convert(T, xf)], tshape...)
    size(xf_arr) == tshape || throw(DimensionMismatch("xf shape $(size(xf_arr)) must match transverse shape $tshape"))

    xnodes = collect(T, CartesianGrids.grid1d(grid, ia))

    if N == 2
        td = _transverse_dims(ia, Val(N))[1]
        snodes_r = CartesianGrids.grid1d(grid, td)
        snodes = collect(T, snodes_r)
        smin, smax = snodes[1], snodes[end]

        xf_vec = vec(xf_arr)
        f_interp = if interp === :linear
            s -> _interp_pwlinear_1d(snodes, xf_vec, s)
        elseif interp === :cubic
            itp = cubic_spline_interpolation(snodes_r, xf_vec, extrapolation_bc=Line())
            s -> itp(s)
        else
            throw(ArgumentError("unsupported interpolation `$interp`; expected :linear or :cubic"))
        end

        @inbounds for I in CartesianIndices(phi)
            s = snodes[I[td]]
            s_eval = periodic ? _wrap_periodic(s, smin, smax) : _clamp_coord(s, smin, smax)
            phi[I] = xnodes[I[ia]] - f_interp(s_eval)
        end
        return phi
    end

    if N == 3
        tdims = _transverse_dims(ia, Val(N))
        ynodes = collect(T, CartesianGrids.grid1d(grid, tdims[1]))
        znodes = collect(T, CartesianGrids.grid1d(grid, tdims[2]))
        ymin, ymax = ynodes[1], ynodes[end]
        zmin, zmax = znodes[1], znodes[end]

        interp === :linear || throw(ArgumentError("3D reconstruction supports interp=:linear"))

        @inbounds for I in CartesianIndices(phi)
            y = ynodes[I[tdims[1]]]
            z = znodes[I[tdims[2]]]
            y_eval = periodic ? _wrap_periodic(y, ymin, ymax) : _clamp_coord(y, ymin, ymax)
            z_eval = periodic ? _wrap_periodic(z, zmin, zmax) : _clamp_coord(z, zmin, zmax)
            phi[I] = xnodes[I[ia]] - _interp_bilinear(ynodes, znodes, xf_arr, y_eval, z_eval)
        end
        return phi
    end

    throw(ArgumentError("phi reconstruction only supports 1D/2D/3D grids"))
end

function phi_from_xf(
    xf,
    grid::CartesianGrids.CartesianGrid{N,T};
    axis::Union{Int,Symbol}=1,
    interp::Symbol=:linear,
    periodic::Bool=false,
) where {N,T}
    phi = Array{T,N}(undef, grid.n...)
    return phi_from_xf!(phi, xf, grid; axis=axis, interp=interp, periodic=periodic)
end
