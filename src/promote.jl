# conversions and promotions from integer and rational types to algebraic number types
import Base.convert
import Base.promote_rule

ExactNumber{T} = Union{T,Rational{T},Complex{T},Complex{Rational{T}}}

# Algebraic number from integer
convert(::Type{AlgebraicNumber{T,F}}, x::ExactNumber{S}) where {S <: Integer,T,F} =
	AlgebraicNumber(promote_type(T, typeof(x))(x), F)

# promotions
promote_rule(::Type{ExactNumber{T}}, ::Type{AlgebraicNumber{S,F}}) where {T <: Integer,S,F} =
	AlgebraicNumber{promote_type(T, S),F}
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
