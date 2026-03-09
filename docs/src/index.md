```@meta
CurrentModule = GlobalHeightFunctions
```

# GlobalHeightFunctions.jl

## Installation

```julia
using Pkg
Pkg.add("GlobalHeightFunctions")
```

## Overview

`GlobalHeightFunctions.jl` is a solver-agnostic toolbox for global height-function interface tracking.

It provides geometry-level operators and does not depend on `PenguinStefan.jl`.

Interfaces must be representable as a single-valued graph along a chosen axis:
- `axis=:x`: `x = xf(y)` in 2D, `x = xf(y,z)` in 3D
- `axis=:y` / `axis=:z`: analogous forms

For multi-valued interfaces along an axis, use a level-set representation instead.

## Quick Example

Start by defining a Cartesian grid and a height function `xf` that represents the interface as a single-valued graph along the `x`-axis.

Then, convert this height function to a signed distance function (SDF) representation using `phi_from_xf`, and convert it back to a height function using `xf_from_sdf`.

```@example
using CartesianGrids
using GlobalHeightFunctions

grid = CartesianGrid((0.0, 0.0), (1.0, 1.0), (65, 49))
ynodes = collect(CartesianGrids.grid1d(grid, 2))
xf = [0.3 + 0.05 * sin(2pi * y) for y in ynodes]

phi = phi_from_xf(xf, grid; axis=:x, interp=:cubic)
xf_hat = xf_from_sdf(phi, grid; axis=:x)
```

## Main Sections

- [API Reference](reference.md)
