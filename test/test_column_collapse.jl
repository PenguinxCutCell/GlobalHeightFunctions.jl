using Test
using CartesianGrids
using GlobalHeightFunctions

@testset "Column collapse and inverse mapping" begin
    grid = CartesianGrid((0.0, 0.0), (1.0, 1.0), (17, 13))
    ynodes = collect(CartesianGrids.grid1d(grid, 2))

    x0 = first(CartesianGrids.grid1d(grid, 1))
    Δy = CartesianGrids.meshsize(grid, 2)
    xf = [0.22 + 0.16 * y for y in ynodes]
    H_expected = (xf .- x0) .* Δy

    V = zeros(Float64, grid.n...)
    @inbounds for j in eachindex(ynodes)
        V[1, j] = H_expected[j]
    end

    Hcol = column_profile(V; axis=1)
    @test maximum(abs.(Hcol .- H_expected)) < 1e-14

    xf_rec = positions_from_heights(Hcol, grid; axis=1)
    @test maximum(abs.(xf_rec .- xf)) < 1e-14
end

@testset "Height update operator" begin
    Hn = [0.1, 0.2, 0.3]
    Hn1 = [0.11, 0.19, 0.31]
    F = [-0.02, 0.0, 0.01]
    rhoL = 2.0

    Hn1_new, R = update_heights(Hn, Hn1, F, rhoL; damping=0.5)

    R_expected = (Hn1 .- Hn) .+ F ./ rhoL
    H_expected = Hn1 .- 0.5 .* R_expected

    @test maximum(abs.(R .- R_expected)) < 1e-14
    @test maximum(abs.(Hn1_new .- H_expected)) < 1e-14
end
