# conversions and promotions from integer and rational types to algebraic number types
import Base.convert
import Base.promote_rule

# Algebraic number from integer
convert(::Type{AlgebraicNumber{T}}, x::S) where {S <: Integer,T} =
	AlgebraicNumber(promote_type(T, S)(x))
convert(::Type{AlgebraicNumber{T}}, x::Rational{S}) where {S <: Integer,T} =
	AlgebraicNumber(Rational{promote_type(T, S)}(x))

convert(::Type{AlgebraicNumber{T}}, x::Complex{S}) where {S <: Integer,T} =
	AlgebraicNumber(Complex{promote_type(T, S)}(x))
convert(::Type{AlgebraicNumber{T}}, x::Complex{Rational{S}}) where {S <: Integer,T} =
	AlgebraicNumber(Complex{Rational{promote_type(T, S)}}(x))

# promotions form Julia types to AlgebraciNumber
promote_rule(::Type{T}, ::Type{AlgebraicNumber{S}}) where {T <: Integer,S} =
	AlgebraicNumber{promote_type(T, S)}
promote_rule(::Type{Rational{T}}, ::Type{AlgebraicNumber{S}}) where {T <: Integer,S} =
	AlgebraicNumber{promote_type(T, S)}

promote_rule(::Type{Complex{T}}, ::Type{AlgebraicNumber{S}}) where {T <: Integer,S} =
	AlgebraicNumber{promote_type(T, S)}
promote_rule(::Type{Complex{Rational{T}}}, ::Type{AlgebraicNumber{S}}) where {T <: Integer,S} =
	AlgebraicNumber{promote_type(T, S)}

# promotion between AlgebraicNumbers
# TODO: it does not work right now
# premote_rule(::Type{AlgebraicNumber{S}}, ::Type{AlgebraicNumber{T}}) where {T,S} =
# 				 AlgebraicNumber{promote_type(T, S)}
# MWE:
# promote(AlgebraicNumber(BigInt(2)), AlgebraicNumber(2))

# conversions back
function convert(::Type{T}, an::AlgebraicNumber) where {T <: Integer}
	if isinteger(an)
		if promote_type(coeffstype(an), T) == T
			return convert(T, -first(moniccoeffs(an)))
		else
			throw(TypeError(:convert, "AlgebraicNumber{$(coeffstype(an)),$(floattype(an))} to $T", coeffstype(an), T))
		end
	else
		throw(InexactError())
	end
end

function convert(::Type{Rational{T}}, an::AlgebraicNumber) where {T <: Integer}
	if isrational(an)
		if promote_type(coeffstype(an), T) == T
			return convert(Rational{T}, -first(moniccoeffs(an)))
		else
			throw(TypeError(:convert, "AlgebraicNumber{$(coeffstype(an)),$(floattype(an))} to Rational{$T}", coeffstype(an), T))
		end
	else
		throw(InexactError())
	end
end

convert(::Type{Rational}, an::AlgebraicNumber) = convert(Rational{coeffstype(an)}, an)

function convert(::Type{Complex{T}}, an::AlgebraicNumber) where {T <: Integer}
	if iscomplexinteger(an)
		if promote_type(coeffstype(an), T) == T
			return convert(Complex{T}, exactdeg2(an))
		else
			throw(TypeError(:convert, "AlgebraicNumber{$(coeffstype(an)),$(floattype(an))} to Complex{$T}", coeffstype(an), T))
		end
		else
		throw(InexactError())
	end
end

function convert(::Type{Complex{Rational{T}}}, an::AlgebraicNumber) where {T <: Integer}
	if iscomplexratioanl(an)
		if promote_type(coeffstype(an), T) == T
			return convert(Complex{Rational{T}}, exactdeg2(an))
		else
			throw(TypeError(:convert, "AlgebraicNumber{$(coeffstype(an)),$(floattype(an))} to Complex{Rational{$T}}", coeffstype(an), T))
		end
		else
		throw(InexactError())
	end
end

function convert(::Type{Complex}, an::AlgebraicNumber)
	if iscomplexinteger(an)
		convert(Complex{coeffstype(an)}, an)
	else
		if iscomplexratioanl(an)
			convert(Complex{Rational{coeffstype(an)}}, an)
		else
		throw(InexactError())
		end
	end
end
