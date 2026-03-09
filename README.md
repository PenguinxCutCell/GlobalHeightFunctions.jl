# GlobalHeightFunctions.jl

`GlobalHeightFunctions.jl` is a solver-agnostic toolbox for **global height-function** interface tracking.

It has no dependency on `PenguinStefan.jl` and only provides geometry-level operators.

## Scope

Interfaces must be representable as a **single-valued graph** along a chosen axis.

- `axis=:x`: `x = xf(y)` in 2D, `x = xf(y,z)` in 3D
- `axis=:y` / `axis=:z`: analogous forms

For multi-valued interfaces along an axis (for example a disk in Cartesian `x=f(y)` form), use a level-set representation instead.

## Provided API

- `xf_from_sdf(phi, grid; axis, method=:zero_crossing)`
- `phi_from_xf!(phi, xf, grid; axis, interp=:linear|:cubic, periodic=false)`
- `phi_from_xf(xf, grid; ...)`
- `column_profile(V; axis)`
- `positions_from_heights(Hcol, grid; axis)`
- `ensure_periodic!(xf)`
- `update_heights(Hn, Hn1, Fcol, rhoL; damping)`

## Quick example

```julia
using CartesianGrids
using GlobalHeightFunctions

grid = CartesianGrid((0.0, 0.0), (1.0, 1.0), (65, 49))
ynodes = collect(CartesianGrids.grid1d(grid, 2))
xf = [0.3 + 0.05*sin(2pi*y) for y in ynodes]

phi = phi_from_xf(xf, grid; axis=:x, interp=:cubic)
xf_hat = xf_from_sdf(phi, grid; axis=:x)
```
