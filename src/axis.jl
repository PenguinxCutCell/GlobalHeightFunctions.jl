@inline _axis_symbol(axis::Symbol) = axis
@inline _axis_symbol(axis::Val{S}) where {S} = S

function axis_to_index(axis, ::Val{N}) where {N}
    if axis isa Integer
        ia = Int(axis)
        1 <= ia <= N || throw(ArgumentError("axis index must be in 1:$N, got $ia"))
        return ia
    end

    sym = _axis_symbol(axis)
    if sym === :x
        return 1
    elseif sym === :y
        N >= 2 || throw(ArgumentError("axis :y requires dimension >= 2"))
        return 2
    elseif sym === :z
        N >= 3 || throw(ArgumentError("axis :z requires dimension >= 3"))
        return 3
    end
    throw(ArgumentError("unsupported axis `$sym`; expected integer, :x, :y, or :z"))
end

@inline axis_to_index(axis, ::Type{A}) where {A<:AbstractArray{<:Any,N}} where {N} = axis_to_index(axis, Val(N))

@inline function _transverse_dims(axis_idx::Int, ::Val{N}) where {N}
    return ntuple(i -> i < axis_idx ? i : i + 1, N - 1)
end

@inline function _line_index(axis_idx::Int, itrans::CartesianIndex{M}, ia::Int, ::Val{N}) where {M,N}
    vals = ntuple(d -> begin
        if d == axis_idx
            ia
        else
            td = d < axis_idx ? d : d - 1
            itrans[td]
        end
    end, N)
    return CartesianIndex(vals)
end
