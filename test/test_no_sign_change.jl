using Test
using CartesianGrids
using GlobalHeightFunctions

@testset "No-sign-change robustness" begin
    grid = CartesianGrid((0.0, 0.0), (1.0, 1.0), (33, 17))

    phi_neg = fill(-1.0, grid.n...)
    xf_neg = xf_from_sdf(phi_neg, grid; axis=:x)
    @test all(isfinite, xf_neg)
    phi_neg_rec = phi_from_xf(xf_neg, grid; axis=:x, interp=:linear)
    @test all(isfinite, phi_neg_rec)

    phi_pos = fill(1.0, grid.n...)
    xf_pos = xf_from_sdf(phi_pos, grid; axis=:x)
    @test all(isfinite, xf_pos)
    phi_pos_rec = phi_from_xf(xf_pos, grid; axis=:x, interp=:linear)
    @test all(isfinite, phi_pos_rec)
end
