
# Define a coeffs iterator for univariate polynomials
struct PolyCoeffs{T <: RingElem}
	poly::T
end

Base.length(p::PolyCoeffs) = max(1, length(p.poly))
Base.eltype(p::PolyCoeffs) = elem_type(base_ring(p.poly))
degree(p::PolyCoeffs) = Nemo.degree(p.poly)

function Base.iterate(p::PolyCoeffs)
	return Nemo.coeff(p.poly, 0), 0
end

function Base.iterate(p::PolyCoeffs, state)
	state += 1
	if degree(p) >= state
		return Nemo.coeff(p.poly, state), state
	else
		return nothing
	end
end

# Now we can add a method to Nemo.coeffs for univariate polynomials
Nemo.coeffs(p::PolyElem) = PolyCoeffs(p)

# Do not restrict to T<: Integer, to use it with fmpz numbers!!
function poly_product_from_coeff(coeffs1::Vector{T}, coeffs2::Vector{T}) where {T}
	R, x = Nemo.PolynomialRing(Nemo.FlintZZ, "x")
	return T.(Nemo.coeffs(R(coeffs1) * R(coeffs2)))
end

poly_product_from_coeff(coeffs1, coeffs2) = poly_product_from_coeff(promote(coeffs1, coeffs2)...)


discriminant_from_coeff(coeffs) = Nemo.discriminant(Nemo.polynomial(Nemo.FlintZZ, coeffs))
Nemo.discriminant(an::AlgebraicNumber) = discriminant_from_coeff(an.minpoly)

function poly_inv_from_coeffs(coeffs, n)
	R, x = Nemo.PowerSeriesRing(Nemo.FlintQQ, n, "x")
	polyinv = Nemo.inv(R(coeffs, length(coeffs), n, 0))
	return [Nemo.coeff(polyinv, i) for i = 0:n - 1]
end

is_integer(x::fmpq) = denominator(x) == 1
Nemo.fmpz(x::fmpq) = is_integer(x) ? numerator(x) : throw(InexactError(:convert, fmpz, x))
Nemo.fmpz(x::Rational) = Base.isinteger(x) ? fmpz(numerator(x)) : throw(InexactError(:convert, fmpz, x))

Base.one(::Type{fmpq}) = fmpq(1)
