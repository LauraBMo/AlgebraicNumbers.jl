
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
