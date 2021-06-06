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

# promotion between AlgebraicNumbers
# TODO: it does not work right now
# premote_rule(::Type{AlgebraicNumber{S,F}}, ::Type{AlgebraicNumber{T,G}}) where {T,F,S,G} =
# 				 AlgebraicNumber{promote_type(T, S),promote_type(F, G)}
# MWE:
# promote(AlgebraicNumber(BigInt(2)), AlgebraicNumber(2))

# conversions back
function convert(::Type{T}, an::AlgebraicNumber) where {T <: Integer}
	if isinteger(an)
		if promote_type(inttype(an), T) == T
			return convert(T, -first(moniccoeffs(an)))
		else
			throw(TypeError(:convert, "AlgebraicNumber{$(inttype(an)),$(floattype(an))} to $T", inttype(an), T))
		end
	else
		throw(InexactError())
	end
end

function convert(::Type{Rational{T}}, an::AlgebraicNumber) where {T <: Integer}
	if isrational(an)
		if promote_type(inttype(an), T) == T
			return convert(Rational{T}, -first(moniccoeffs(an)))
		else
			throw(TypeError(:convert, "AlgebraicNumber{$(inttype(an)),$(floattype(an))} to Rational{$T}", inttype(an), T))
		end
	else
		throw(InexactError())
	end
end

convert(::Type{Rational}, an::AlgebraicNumber) = convert(Rational{inttype(an)}, an)

function convert(::Type{Complex{T}}, an::AlgebraicNumber) where {T <: Integer}
	if iscomplexinteger(an)
		if promote_type(inttype(an), T) == T
			return convert(Complex{T}, exactdeg2(an))
		else
			throw(TypeError(:convert, "AlgebraicNumber{$(inttype(an)),$(floattype(an))} to Complex{$T}", inttype(an), T))
		end
		else
		throw(InexactError())
	end
end

function convert(::Type{Complex{Rational{T}}}, an::AlgebraicNumber) where {T <: Integer}
	if iscomplexratioanl(an)
		if promote_type(inttype(an), T) == T
			return convert(Complex{Rational{T}}, exactdeg2(an))
		else
			throw(TypeError(:convert, "AlgebraicNumber{$(inttype(an)),$(floattype(an))} to Complex{Rational{$T}}", inttype(an), T))
		end
		else
		throw(InexactError())
	end
end

function convert(::Type{Complex}, an::AlgebraicNumber)
	if iscomplexinteger(an)
		convert(Complex{inttype(an)}, an)
	else
		if iscomplexratioanl(an)
			convert(Complex{Rational{inttype(an)}}, an)
		else
		throw(InexactError())
		end
	end
end
