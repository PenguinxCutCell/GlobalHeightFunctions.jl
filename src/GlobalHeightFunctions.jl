module GlobalHeightFunctions

using LinearAlgebra
using StaticArrays

using CartesianGrids
using Interpolations

export axis_to_index
export spatial_shape_from_grid, reshape_V, column_profile, positions_from_heights, ensure_periodic!
export xf_from_sdf
export phi_from_xf!, phi_from_xf
export update_heights

include("axis.jl")
include("height_tracking.jl")
include("sdf_to_height.jl")
include("interp.jl")
include("height_to_sdf.jl")
include("update.jl")

end # module
