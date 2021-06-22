# _Fast_ composed sums and composed products of polynomials,
# using the algorithm described in:
# "Fast computation of special resultants"
# by Bostan, Flajolet, Salvy, and Schost

# compute newton power series of polynomial given with coefficients coeff,
# in base field R,x.
# See fig.1 in reference
function poly_to_Newton(coeffs::Vector{T}, N, R, x) where {T}
	# first, make monic.
	mcoeffs = moniccoeffs(coeffs)
	D = degree(mcoeffs)

	# initialize power series polynomials
	A = R(reverse(derivative(mcoeffs)))
	B = R(reverse(mcoeffs))

	b0 = R(poly_inv_from_coeffs(reverse(mcoeffs), D))
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
	return Nemo.coeffs(r)
end

# This algorithm is based on the Leverrier-Faddeev algorithm
# see: http://math.stackexchange.com/questions/405822/what-is-the-fastest-way-to-find-the-characteristic-polynomial-of-a-matrix
# using LinearAlgebra: dot
# Old implematation
# function from_newton{T}(tr::Vector{T})
# 	c = T[]
# 	for k = 1 : length(tr)-1
# 		push!(c, -dot(tr[2:(k+1)], vcat(reverse(c),1))/k)
# 	end
# 	return vcat(reverse(c),1)
# end

function Newton_to_poly(N::Vector{T}, D=length(N)) where {T}
	# special case
	n = length(N)
	out = T[0,1]
	if N != [1]
		c = zeros(T, max(D, n)) # The first D~n entries are zero.
		c[end] = one(T)
		for k = 1:n - 1
			next_c = -sum(N[2:(k + 1)] .* c[(end - k + 1):end]) // k
			c[end - k] = next_c
		end
		out = c
	end
	return out
end
