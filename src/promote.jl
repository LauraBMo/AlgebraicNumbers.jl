# conversions and promotions from integer and rational types to algebraic number types
import Base.convert
import Base.promote_rule

# Algebraic number from integer
convert(::Type{AlgebraicNumber{T,F}}, x::S) where {S <: Integer,T,F} =
	AlgebraicNumber(promote_type(T, S)(x), F)
convert(::Type{AlgebraicNumber{T,F}}, x::Rational{S}) where {S <: Integer,T,F} =
	AlgebraicNumber(Rational{promote_type(T, S)}(x), F)
convert(::Type{AlgebraicNumber{T,F}}, x::Complex{S}) where {S <: Integer,T,F} =
	AlgebraicNumber(Complex{promote_type(T, S)}(x), F)
convert(::Type{AlgebraicNumber{T,F}}, x::Complex{Rational{S}}) where {S <: Integer,T,F} =
	AlgebraicNumber(Complex{Rational{promote_type(T, S)}}(x), F)

# promotions form Julia types to AlgebraciNumber
promote_rule(::Type{T}, ::Type{AlgebraicNumber{S,F}}) where {T <: Integer,S,F} =
	AlgebraicNumber{promote_type(T, S),F}
promote_rule(::Type{Rational{T}}, ::Type{AlgebraicNumber{S,F}}) where {T <: Integer,S,F} =
	AlgebraicNumber{promote_type(T, S),F}
promote_rule(::Type{Complex{T}}, ::Type{AlgebraicNumber{S,F}}) where {T <: Integer,S,F} =
	AlgebraicNumber{promote_type(T, S),F}
promote_rule(::Type{Complex{Rational{T}}}, ::Type{AlgebraicNumber{S,F}}) where {T <: Integer,S,F} =
	AlgebraicNumber{promote_type(T, S),F}
premote_rule(::Type{AlgebraicNumber{T,F}},::Type{AlgebraicNumber{S,G}}) where {T,S <: Integer,F,G <: AbstractFloat} =
	AlgebraicNumber{promote_type(T, S),promote_type(F, G)}

# promotion between AlgebraicNumbers
premote_rule(::Type{AlgebraicNumber{T,F}}, ::Type{AlgebraicNumber{S,G}}) where {T,S,F,G} =
	AlgebraicNumber{promote_type(T, S),promote_type(F, G)}

# conversions back
function convert(::Type{T}, an::AlgebraicNumber) where {T <: Integer}
	if isinteger(an)
		return convert(T, -prod(an.coeffs))
	else
		throw(InexactError())
	end
end

function convert(::Type{Rational{T}}, an::AlgebraicNumber) where {T <: Integer}
	if isrational(an)
		return convert(Rational{T}, -first(moniccoeffs(an)))
	else
		throw(InexactError())
	end
end
