

# derivative of polynomial
degree(coeffs::Vector) = length(coeffs) - 1
degree(an::AlgebraicNumber) = degree(an.coeffs)
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
