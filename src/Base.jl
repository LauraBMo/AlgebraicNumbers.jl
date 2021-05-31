

# derivative of polynomial
degree(coeffs::Vector) = length(coeffs) - 1
derivative(coeffs::Vector) = coeffs[2:end] .* (1:degree(coeffs))

# get monicoeffs for minpoly from an an
moniccoeffs(coeffs::Vector) = coeffs // coeffs[end]
moniccoeffs(an::AlgebraicNumber) = moniccoeffs(an.coeff)[begin:(end - 1)]

# if x is root of coeff, then -x is root of minus_minpoly(coeff)
minus_minpoly(coeff::Vector) =  [isodd(i) ? -coeff[i] : coeff[i] for i = 1:length(coeff)]


# given an array of Rationals computes the minimal integer multiple
convert_intcoeffs(::Type{T}, v) where {T <: Integer} = T.(lcm(denominator.(v)) * v)
convert_intcoeffs(v) = convert_intcoeffs(BigInt, v)
