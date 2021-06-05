
import Base: hash, show

# An algebraic number,
# consisting of the minimal polynomial of the number,
# an arbitrary-precision approximation of the number,
# and prec which specifies the minimal distance between roots of p
# TODO: apprx has to be complex.
"""
        AlgebraicNumber{T}

An algebraic number consisting of
   - the minimal polynomial of the number: Vector of integers.
   - Complex float point aproximation
   - minimal distance between two roots of minimal polynomial
"""
struct AlgebraicNumber{T <: Integer,F <: AbstractFloat} <: Number
    coeffs::Vector{T}
    apprx::Complex{F}
    prec::F
end

Base.hash(an::AlgebraicNumber, h::UInt) = hash((an.coeffs, an.apprx), h)

# TODO: only show up to precision
function Base.show(io::IO, an::AlgebraicNumber)
    print(io, "≈")
    # ndigits = max(10, round(Int,ceil(convert(Float64,log(an.prec)/log(10)))))
    show(io, convert(Complex{Float64}, an.apprx))
    # print(io,"...")
end
