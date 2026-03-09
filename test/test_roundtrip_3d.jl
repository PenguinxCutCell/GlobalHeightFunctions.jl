using Test
using CartesianGrids
using GlobalHeightFunctions

function _build_phi_from_xf_3d(grid, xf_fun)
    xnodes = collect(CartesianGrids.grid1d(grid, 1))
    ynodes = collect(CartesianGrids.grid1d(grid, 2))
    znodes = collect(CartesianGrids.grid1d(grid, 3))
    phi = Array{Float64}(undef, grid.n...)
    @inbounds for k in eachindex(znodes), j in eachindex(ynodes), i in eachindex(xnodes)
        phi[i, j, k] = xnodes[i] - xf_fun(ynodes[j], znodes[k])
    end
    return phi
end

@testset "SDF <-> Global HF roundtrip (3D smoke)" begin
    xf_fun = (y, z) -> 0.42 + 0.03 * sin(2pi * y) * sin(2pi * z)

    errs = Float64[]
    for n in (17, 33)
        grid = CartesianGrid((0.0, 0.0, 0.0), (1.0, 1.0, 1.0), (n, n, n))
        phi = _build_phi_from_xf_3d(grid, xf_fun)
        xf = xf_from_sdf(phi, grid; axis=:x)
        phi_hat = phi_from_xf(xf, grid; axis=:x, interp=:linear)
        err = maximum(abs.(phi_hat .- phi))
        push!(errs, err)
        @test isfinite(err)
    end

    @test errs[2] <= errs[1] + 1e-12
end
