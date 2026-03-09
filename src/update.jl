function update_heights(
    Hn,
    Hn1,
    Fcol,
    rhoL;
    damping::Real=0.5,
)
    size(Hn) == size(Hn1) || throw(DimensionMismatch("Hn and Hn1 must have identical shapes"))
    size(Hn) == size(Fcol) || throw(DimensionMismatch("Fcol shape must match Hn"))

    T = promote_type(eltype(Hn), eltype(Hn1), eltype(Fcol), typeof(rhoL), typeof(damping))
    ρL = convert(T, rhoL)
    ω = convert(T, damping)

    Rcol = (Hn1 .- Hn) .+ Fcol ./ ρL
    Hn1_new = Hn1 .- ω .* Rcol
    return Hn1_new, Rcol
end
