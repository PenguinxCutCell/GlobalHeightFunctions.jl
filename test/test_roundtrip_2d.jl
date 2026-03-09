using Test
using CartesianGrids
using GlobalHeightFunctions

@testset "SDF <-> Global HF roundtrip" begin
    grid = CartesianGrid((0.0, 0.0), (1.0, 1.0), (65, 49))
    xnodes = collect(CartesianGrids.grid1d(grid, 1))
    ynodes = collect(CartesianGrids.grid1d(grid, 2))

    xf_exact = [0.37 + 0.08 * sin(2pi * y) for y in ynodes]
    phi = Array{Float64}(undef, grid.n...)
    @inbounds for j in eachindex(ynodes), i in eachindex(xnodes)
        phi[i, j] = xnodes[i] - xf_exact[j]
    end

    xf_hat = xf_from_sdf(phi, grid; axis=:x)
    ensure_periodic!(xf_hat)
    @test maximum(abs.(xf_hat .- xf_exact)) < 5e-13
    @test xf_hat[end] == xf_hat[1]

    phi_lin = phi_from_xf(xf_hat, grid; axis=:x, interp=:linear)
    phi_cub = phi_from_xf(xf_hat, grid; axis=:x, interp=:cubic)

    @test maximum(abs.(phi_lin .- phi)) < 5e-13
    @test maximum(abs.(phi_cub .- phi)) < 5e-13
end
