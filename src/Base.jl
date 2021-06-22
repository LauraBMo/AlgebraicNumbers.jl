
# Degree of a polynomial as a vector of coefficients
degree(coeffs::Vector) = length(coeffs) - 1
degree(an::AlgebraicNumber) = degree(an.minpoly)

# derivative of polynomial
derivative(coeffs::Vector) = coeffs[2:end] .* (1:degree(coeffs))

# get monicoeffs for minpoly from an an
moniccoeffs(coeffs::Vector) = coeffs .// coeffs[end]
moniccoeffs(an::AlgebraicNumber) = moniccoeffs(an.minpoly)
ismonic(an::AlgebraicNumber) = abs(an.minpoly[end]) == 1

# if x is root of coeff, then -x is root of minus_minpoly(coeff)
minus_minpoly(coeffs::Vector) =  [isodd(i) ? -coeffs[i] : coeffs[i] for i = 1:length(coeffs)]


vec_to_poly(coeffs::Vector, var) = dot(coeffs, var.^(0:degree(coeffs)))
confirm_algnumber(an::AlgebraicNumber) = vec_to_poly(big.(an.minpoly), big.(an.approx))

isrational(an::AlgebraicNumber) = degree(an) == 1
isinteger(an::AlgebraicNumber) = isrational(an) && ismonic(an)

function isdiscriminantsquare(an::AlgebraicNumber)
	Δ = Nemo.discriminant(an)
	return issquare(abs(Δ))
end

iscomplexratioanl(an::AlgebraicNumber) = degree(an) == 2 && isdiscriminantsquare(an::AlgebraicNumber)
iscomplexinteger(an::AlgebraicNumber) = iscomplexratioanl(an) && ismonic(an)

sqrtdiscriminant(an::AlgebraicNumber) = isqrt(abs(Nemo.discriminant(an)))

# ±(a,b) = (a + b, a - b)
# rootsdeg2(an::AlgebraicNumber) = (-an.coeffs[2] ± sqrtdiscriminant(an) * Complex{inttype(an)}(0, 1)) .// 2

function exactdeg2(::Type{T}, an::AlgebraicNumber) where {T}
	coeffs = T.(an.minpoly)
	Δ = T(sqrtdiscriminant(an))
	if imag(an.approx) > 0
		num = (-coeffs[2] + Δ * Complex{T}(0, 1)) // (2 * coeffs[3])
	else
		num = (-coeffs[2] - Δ * Complex{T}(0, 1)) // (2 * coeffs[3])
	end
	return num
end

exactdeg2(an) = exactdeg2(BigInt, an)

"""
    integertype(::Type{S})

Returns the integer type for Julia exact numbers.
That is, `T` for `S` in `Union{T, Rational{T}, Complex{T}, Complex{Rational{T}}}`

There are also the methods

    integertype(x::Number) = integertype(typeof(x))
    integertype(x::Array) = integertype(eltype(x))

"""
function integertype end

integertype(::Type{T}) where {T <: Integer} = T
integertype(::Type{Rational{T}}) where {T <: Integer} = T
integertype(::Type{Complex{T}}) where {T <: Integer} = T
integertype(::Type{Complex{Rational{T}}}) where {T <: Integer} = T

integertype(x::Number) = integertype(typeof(x))
integertype(x::Array) = integertype(eltype(x))
