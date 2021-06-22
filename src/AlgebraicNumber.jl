
import Base: hash, show, convert, promote_rule

"""
        AlgebraicNumber

Structure for algebraic numbers. It consists of
   - the minimal polynomial of the number, a vector of fmpz,
   - complex ball containing a certified root of minpoly.


     AlgebraicNumber(x)

Converts `x` to `AlgebraicNumber`, where `x` is of some Julia exact number type,
`Union{S, Rational{S}, Complex{S}, Complex{Rational{S}}} where {S <: Integer}`.
"""
struct AlgebraicNumber <: Number
    minpoly::Vector{fmpz}
    approx::Acb
    # Arblib uses RefVal DEFAULT_PRECISION[] to set new Arb/Acb....
end

getminpoly(an::AlgebraicNumber) = an.minpoly
getapprox(an::AlgebraicNumber) = an.approx

Base.hash(an::AlgebraicNumber, h::UInt) = hash((an.minpoly, an.apprx), h)
Base.show(io::IO, an::AlgebraicNumber) = print(io, "≃" * Arblib.string_nice(an.approx, 20))

const IntOrRat = Union{Integer,Rational}

# Create Algebraic Numbers from Julia exact numbers.
AlgebraicNumber(x::T) where {T <: IntOrRat} = AlgebraicNumber([-numerator(x), denominator(x)], x)

function AlgebraicNumber(x::Complex{T}) where {T <: Integer}
    return isreal(x) ? AlgebraicNumber(real(x)) :
        AlgebraicNumber([imag(x)^2 + real(x)^2, -2 * real(x), one(T)], x)
end

function AlgebraicNumber(x::Complex{Rational{T}}) where {T <: Integer}
    return isreal(x) ? AlgebraicNumber(real(x)) :
        AlgebraicNumber(RatVec([imag(x)^2 + real(x)^2, -2 * real(x), one(T)]), x)
end

# conversions and promotions from integer and rational types to algebraic number types

# Algebraic number from integer
# using duck catch
convert(::Type{AlgebraicNumber}, x::T) where {T <: Number} = AlgebraicNumber(x)

# promotions form Julia types to AlgebraciNumber
promote_rule(::Type{T}, ::Type{AlgebraicNumber}) where {T <: Number} = AlgebraicNumber

# TODO: it does not work right now
# premote_rule(::Type{AlgebraicNumber{S}}, ::Type{AlgebraicNumber{T}}) where {T,S} =
# 				 AlgebraicNumber{promote_type(T, S)}
# MWE:
# promote(AlgebraicNumber(BigInt(2)), AlgebraicNumber(2))

# conversions back
function convert(::Type{T}, an::AlgebraicNumber) where {T <: IntOrRat}
    return isrational(an) ?
        -first(moniccoeffs(integertype(T).(an.minpoly))) :
        throw(InexactError(:convert, T, an))
end

convert(::Type{Rational}, an::AlgebraicNumber) = convert(Rational{BigInt}, an)

function convert(::Type{Complex{T}}, an::AlgebraicNumber) where {T <: IntOrRat}
    if isrational(an)
        return convert(Complex{T}, convert(T, an))
    else
        return iscomplexratioanl(an) ?
            # If iscomplexinteger, eactdeg2 will return numbers convertibles to int.
            exactdeg2(integertype(T), an) :
            throw(InexactError(:convert, T, an))
        end
end

function convert(::Type{Complex}, an::AlgebraicNumber)
    if isrational(an)
        if isinteger(an)
            return convert(Complex{BigInt}, convert(BigInt, an))
        else
            return convert(Complex{Rational{BigInt}}, convert(Rational, an))
        end
    else
        if iscomplexinteger(an)
            return convert(Complex{BigInt}, an)
        else
            if iscomplexratioanl(an)
                return convert(Complex{Rational{BigInt}}, an)
            else
                throw(InexactError(:convert, T, an))
            end
        end
    end
end
