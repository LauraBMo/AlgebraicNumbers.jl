

# # Hadamard (element-wise) product of two polynomials

"""
    factorials(n)

List of the n first factorial numbers, the coefficients of exponential series.
"""
factorials(n) = factorial.(Nemo.FlintZZ.(0:n))
coeffsexp(n) = 1 .// factorials(n)
coeffsinvexp(n) = factorials(n)

Hadamard(coeffs1, coeffs2) = prod.(zip(coeffs1, coeffs2))

# composed product of two polynomials, given as coeffs p and q
function composed_product(coeffsp::Vector{T}, coeffsq::Vector{T}) where {T}
	# compute newton series
	n = degree(coeffsp) * degree(coeffsq) + 1
	R, x = Nemo.PolynomialRing(Nemo.FlintQQ, "x")
	# eats coeffsp::Vector{fmpz}; returns iterator coeffs of poly in R
	coeffsA = poly_to_Newton(coeffsp, n, R, x)
	coeffsB = poly_to_Newton(coeffsq, n, R, x)

	# Hadamard product of the newton series and convert back to polynomial
	coeffsAB = Hadamard(coeffsA, coeffsB)
	coeffspq = Newton_to_poly(coeffsAB, n + 1)

	# convert to integer and return
	return convert(Vector{T}, RatVec(coeffspq))
end

composed_product(p,q) = composed_product(promote(p, q)...)

# composed sum of two polynomials, given as coeffs p and q
function composed_sum(coeffsp::Vector{T}, coeffsq::Vector{T}) where {T}
	# compute newton series
	n = degree(coeffsp) * degree(coeffsq) + 1
	R, x = Nemo.PolynomialRing(Nemo.FlintQQ, "x")
	coeffsA = poly_to_Newton(coeffsp, n, R, x)
	coeffsB = poly_to_Newton(coeffsq, n, R, x)

	# Multiply newton series
	m = Nemo.coeffs(mullow(
		R(Hadamard(coeffsA, coeffsexp(degree(coeffsA)))),
		R(Hadamard(coeffsB, coeffsexp(degree(coeffsB)))),
		n + 1
	))
	# Convert back to polynomial
	pq = Newton_to_poly(Hadamard(m, coeffsinvexp(degree(m))), n + 1)
	# convert to integer and return
	return convert(Vector{T}, RatVec(pq))
end
