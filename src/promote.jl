# conversions and promotions from integer and rational types to algebraic number types
import Base.convert
import Base.promote_rule

# Algebraic number from integer
convert(::Type{T}, x::Number) where {T <: AlgebraicNumber} = AlgebraicNumber(x)

# promotions
promote_rule(::Type{T}, ::Type{AlgebraicNumber{S,F}}) where {T <: Number,S,F} = AlgebraicNumber{S,F}

isrational(an::AlgebraicNumber) = degree(an.coeff) == 1
isinteger(an::AlgebraicNumber) = isrational(an) && abs(an.coeff[2]) == 1

# conversions back
function convert(::Type{T}, an::AlgebraicNumber) where {T <: Integer}
	if isinteger(an)
		return convert(T, -prod(an.coeff))
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
