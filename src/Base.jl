

# Degree of a polynomial as a vector of coefficients
degree(coeffs::Vector) = length(coeffs) - 1
degree(an::AlgebraicNumber) = degree(an.coeffs)

# derivative of polynomial
derivative(coeffs::Vector) = coeffs[2:end] .* (1:degree(coeffs))

# get monicoeffs for minpoly from an an
moniccoeffs(coeffs::Vector) = coeffs // coeffs[end]
moniccoeffs(an::AlgebraicNumber) = moniccoeffs(an.coeffs)[begin:(end - 1)]
ismonic(an::AlgebraicNumber) = abs(an.coeffs[end]) == 1

# if x is root of coeff, then -x is root of minus_minpoly(coeff)
minus_minpoly(coeffs::Vector) =  [isodd(i) ? -coeffs[i] : coeffs[i] for i = 1:length(coeffs)]


# given an array of Rationals computes the minimal integer multiple
convert_intcoeffs(::Type{T}, v) where {T <: Integer} = T.(lcm(denominator.(v)) * v)
convert_intcoeffs(v) = convert_intcoeffs(BigInt, v)

confirm_algnumber(an) = sum(an.coeffs .* (an.apprx.^(0:degree(an.coeffs))))

function iswelldefined(mindist::Real, prec::Real)
	mindist < prec || throw("Error!!! apprx is a distance bigger than prec to any root.")
end

iswelldefined(coeff::Vector{T}, apprx, prec) where {T <: Integer} = iswelldefined(distance(coeff, apprx), prec)
iswelldefined(an::AlgebraicNumber) = iswelldefined(an.coeffs, an.apprx, an.prec)

isrational(an::AlgebraicNumber) = degree(an) == 1
isinteger(an::AlgebraicNumber) = isrational(an) && ismonic(an)

function isdiscriminantsquare(an::AlgebraicNumber)
	Δ = Nemo.discriminant(an)
	return issquare(Δ) || issquare(-Δ)
end

iscomplexratioanl(an::AlgebraicNumber) = degree(an) == 2 && isdiscriminantsquare(an::AlgebraicNumber)
iscomplexinteger(an::AlgebraicNumber) = isratcomplex(an) && ismonic(an)

function rootscomplexrational(an::AlgebraicNumber)
	mcoeffs = moniccoeffs(an)
end

sqrtdiscriminant(an::AlgebraicNumber) = inttype(an)(isqrt(-Nemo.discriminant(an)))

# ±(a,b) = (a + b, a - b)
# rootsdeg2(an::AlgebraicNumber) = (-an.coeffs[2] ± sqrtdiscriminant(an) * Complex{inttype(an)}(0, 1)) .// 2

function exactdeg2(an::AlgebraicNumber)
	if im(an.apprx) > 0
		num = (-an.coeffs[2] + sqrtdiscriminant(an) * Complex{inttype(an)}(0, 1)) .// 2
	else
		num = (-an.coeffs[2] - sqrtdiscriminant(an) * Complex{inttype(an)}(0, 1)) .// 2
	end
	return num
end

# AlgebraicNumbers.solvedeg2(x)
