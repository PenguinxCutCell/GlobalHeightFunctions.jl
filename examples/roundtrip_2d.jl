using CartesianGrids
using GlobalHeightFunctions

grid = CartesianGrid((0.0, 0.0), (1.0, 1.0), (65, 49))
ynodes = collect(CartesianGrids.grid1d(grid, 2))
xf = [0.3 + 0.05 * sin(2pi * y) for y in ynodes]

phi = phi_from_xf(xf, grid; axis=:x, interp=:cubic)
xf_hat = xf_from_sdf(phi, grid; axis=:x)

println("max |xf_hat - xf| = ", maximum(abs.(xf_hat .- xf)))
