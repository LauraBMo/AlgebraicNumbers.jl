
"""
    RatVec{T}

Vectors of eltype `T`.

Used to avoid type piracy in converting them to vectors
of some integer of type multiplying by the lcm of the denominators.
"""
struct RatVec{T} <: AbstractVector{T}
	vec::Vector{T}
end

getvector(v::RatVec) = v.vec

Base.length(v::RatVec) = Base.length(getvector(v))
Base.size(v::RatVec) = Base.size(getvector(v))

Base.eltype(v::RatVec{T}) where T = T

Base.IndexStyle(::Type{<:RatVec}) = Base.IndexLinear()
Base.getindex(v::RatVec, i::Int) = Base.getindex(getvector(v), i)

Base.convert(::Type{Vector{T}}, v::RatVec{S}) where {T,S} =
    T.(lcm(denominator.(getvector(v))) .* getvector(v))
