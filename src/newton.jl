# _Fast_ composed sums and composed products of polynomials,
# using the algorithm described in:
# "Fast computation of special resultants"
# by Bostan, Flajolet, Salvy, and Schost

function polyinv(coeffs::Vector, n)
	R, x = Nemo.PowerSeriesRing(Nemo.FlintQQ, n, "x")
	a = R(coeffs, length(coeffs), n, 0)
	ai = inv(a)
	return [coeff(ai, i) for i = 0:n - 1]
end

# compute newton power series of polynomial given with coefficients coeff,
# in base field R,x.
# See fig.1 in reference
function poly_to_newton(coeffs::Vector{T}, N, R, x) where {T <: Integer}
	# first, make monic.
	mcoeffs = moniccoeffs(coeffs)
	D = degree(mcoeffs)

	# initialize power series polynomials
	A = R(reverse(derivative(mcoeffs)))
	B = R(reverse(mcoeffs))

	b0 = R(polyinv(reverse(mcoeffs), D))
	c  = mullow(A, b0, D)

	r = R()
	x_power = R(1)
	x_d = x^D

	step(C) = -mullow(shift_right(B * C, D), b0, D)
	for j = 0:Int(floor(N / D))
		r += c * x_power
		x_power *= x_d
		c = step(c)
	end
	return Rational.(Nemo.coeffs(r))
end

# to_array(p) = Rational{BigInt}[Rational(coeff(p, i)) for i = 0:Nemo.degree(p)]

# tr: traces i.e. newton series
# This algorithm is based on the Leverrier-Faddeev algorithm
# see: http://math.stackexchange.com/questions/405822/what-is-the-fastest-way-to-find-the-characteristic-polynomial-of-a-matrix
# using LinearAlgebra: dot
# function from_newton{T}(tr::Vector{T})
# 	c = T[]
# 	for k = 1 : length(tr)-1
# 		push!(c, -dot(tr[2:(k+1)], vcat(reverse(c),1))/k)
# 	end
# 	return vcat(reverse(c),1)
# end 

function newton_to_poly(N::Vector{T}, D=length(N)) where {T <: Number}
	# special case
	n = length(N)
	out = T[0,1]
	if N != [1]
		c = Vector{T}(undef, n)
		c[end] = one(T)
		for k = 1:n - 1
			next_c = -sum(N[2:(k + 1)] .* c[(end - k + 1):end]) // k
			c[end - k] = next_c
		end
		out = c
	end
	while length(out) < D
		pushfirst!(out, zero(T))
	end
	return out
end

# Hadamard (element-wise) product of two polynomials
function Hadamard(p, q)
	return Rational{BigInt}[cp * cq for (cp, cq) in zip(Nemo.coeffs(p), Nemo.coeffs(q))]
end

Hadamard_exp(v) = v.*(1 .// factorial.(big.(0:degree(v))))
# Hadamard_invexp(v) =v.*(factorial.(big.(0:degree(v))))

# Hadamard_exp(p) = [Rational(c) // factorial(big(i - 1)) for (i, c) in enumerate(Nemo.coeffs(p))]
Hadamard_invexp(p) = [Rational(c) * factorial(big(i - 1)) for (i, c) in enumerate(Nemo.coeffs(p))]

# composed product of two polynomials, given as coeffs p and q
function composed_product(p::Vector{T}, q::Vector{T}) where {T <: Integer}
	# compute newton series
	n = degree(p) * degree(q) + 1
	R, x = Nemo.PolynomialRing(Nemo.FlintQQ, "x")
	a = poly_to_newton(p, n, R, x)
	b = poly_to_newton(q, n, R, x)

	# multiply newton series and invert
	pq = newton_to_poly(prod.(zip(a,b)), n)

	# convert to integer and return
	return convert_intcoeffs(T, pq)
end

composed_product(p,q) = composed_product(promote(p, q)...)

# composed sum of two polynomials, given as coeffs p and q
function composed_sum(p::Vector{T}, q::Vector{T}) where {T <: Integer}
	# compute newton series
	n = degree(p) * degree(q) + 1
	R, x = Nemo.PolynomialRing(Nemo.FlintQQ, "x")
	a = poly_to_newton(p, n, R, x)
	b = poly_to_newton(q, n, R, x)

	# multiply newton series and invert
	m = mullow(R(Hadamard_exp(a)), R(Hadamard_exp(b)), n + 1)
	pq = newton_to_poly(Hadamard_invexp(m), n)
	# convert to integer and return
	return convert_intcoeffs(T, pq)
end

composed_sum(p,q) = composed_sum(promote(p, q)...)
