# _Fast_ composed sums and composed products of polynomials,
# using the algorithm described in:
# "Fast computation of special resultants"
# by Bostan, Flajolet, Salvy, and Schost

# derivative of polynomial
degree(v::Vector) = length(v) - 1
derivative(c::Vector) = c[2:end] .* (1:degree(c))

function polyinv(coeffs::Vector, n)
	R, x = Nemo.PowerSeriesRing(Nemo.FlintQQ, n, "x")
	a = R(coeffs, length(coeffs), n, 0)
	ai = inv(a)
	return [coeff(ai, i) for i = 0:n - 1]
end

# compute newton power series of polynomial given with coefficients coeff,
# in base field R,x.
# See fig.1 in reference
function poly_to_newton(coeffs::Vector{T}, n, R, x) where {T <: Integer}
	# first, make monic.
	coeffs = coeffs // coeffs[end]

	d = degree(coeffs)
	a_cfs = reverse(derivative(coeffs))
	b_cfs = reverse(coeffs)

	# initialize power series polynomials
	a = R(a_cfs)
	b = R(b_cfs)
	b0 = R(polyinv(b_cfs, n))

	c  = truncate(a * b0, d)

	r = R()
	x_power = R(1)
	x_d = x^d

	l = Int(floor(n / d))
	for j = 0:l
		r += c * x_power
		x_power *= x_d
		c = -mullow(shift_right(b * c, d), b0, d)
	end
	return r
end

# to_array(p) = Rational{BigInt}[Rational(coeff(p, i)) for i = 0:Nemo.degree(p)]

# tr: traces i.e. newton series
# This algorithm is based on the Leverrier-Faddeev algorithm
# see: http://math.stackexchange.com/questions/405822/what-is-the-fastest-way-to-find-the-characteristic-polynomial-of-a-matrix
function newton_to_poly(tr::Vector{T}) where {T <: Number}
	# special case
	out = T[0,1]
	if tr != [1]
		n = length(tr)
		c = Vector{T}(undef, n)
		c[end] = one(T)
		for k = 1:n - 1
			next_c = -sum(tr[2:(k + 1)] .* c[(end - k + 1):end]) // k
			c[end - k] = next_c
		end
		out = c
	end
	return out
end

# Hadamard (element-wise) product of two polynomials
function Hadamard(p, q)
	n = min(Nemo.degree(p), Nemo.degree(q))
	return Rational{BigInt}[Nemo.coeff(p, i) * Nemo.coeff(q, i) for i = 0:n]
end

# composed product of two polynomials, given as coeffs p and q
function composed_product(p::Vector{T}, q::Vector{T}) where {T <: Integer}
	# compute newton series
	n = degree(p) * degree(q) + 1
	R, x = Nemo.PolynomialRing(Nemo.FlintQQ, "x")
	a = poly_to_newton(p, n, R, x)
	b = poly_to_newton(q, n, R, x)

	# multiply newton series and invert
	pq = newton_to_poly(Hadamard(a, b))

	# convert to integer and return
	return convert_intcoeffs(T, pq)
end

# composed sum of two polynomials, given as coeffs p and q
function composed_sum(p::Vector{T}, q::Vector{T}) where {T <: Integer}
	# compute newton series
	n = degree(p) * degree(q) + 1
	R, x = Nemo.PolynomialRing(Nemo.FlintQQ, "x")
	a = poly_to_newton(p, n, R, x)
	b = poly_to_newton(q, n, R, x)

	# exp series
	ee  = R([Nemo.FlintQQ(1 // factorial(BigInt(i))) for i = 0:n])
	eei = R([Nemo.FlintQQ(factorial(BigInt(i))) for i = 0:n])

	# multiply newton series and invert
	m = mullow(R(Hadamard(a, ee)), R(Hadamard(b, ee)), n + 1)
	pq = newton_to_poly(Hadamard(m, eei))

	# convert to integer and return
	return convert_intcoeffs(T, pq)
end

convert_intcoeffs(::Type{T}, v) where {T <: Integer} = T.(lcm(denominator.(v)) * v)
