
# Define a coeffs iterator for univariate polynomials
struct PolyCoeffs{T <: RingElem}
	poly::T
end

Base.length(p::PolyCoeffs) = max(1, length(p.poly))
Base.eltype(p::PolyCoeffs) = elem_type(base_ring(p.poly))

function Base.iterate(x::PolyCoeffs)
	return Nemo.coeff(x.poly, 0), 0
end

function Base.iterate(x::PolyCoeffs, state)
	state += 1
	if Nemo.degree(x.poly) >= state
		return Nemo.coeff(x.poly, state), state
	else
		return nothing
	end
end

# Now we can add a method to Nemo.coeffs for univariate polynomials
Nemo.coeffs(p::PolyElem) = PolyCoeffs(p)

function poly_product_from_coeff(coeffs1, coeffs2)
	R, x = PolynomialRing(Nemo.FlintZZ, "x")
	return get_coeffs(R(coeffs1) * R(coeffs2), promote_type(eltype(coeffs1), eltype(coeffs2)))
end

discriminant_from_coeff(coeffs) = Nemo.discriminant(Nemo.polynomial(Nemo.FlintZZ, coeffs))
Nemo.discriminant(an::AlgebraicNumber) = discriminant_from_coeff(an.coeffs)
