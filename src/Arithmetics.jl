# Exact representation of algebraic numbers
# (Numbers that are roots of polynomials with integer coefficients)
# And also arithmetic on algebraic numbers,
# including +, -, *, /, and radicals.

import Base: zero, one, abs, conj, real, imag
import Base: +, -, *, /, inv, ==, sqrt, cbrt

function ==(an1::AlgebraicNumber, an2::AlgebraicNumber)
	return moniccoeffs(an1) == moniccoeffs(an2) &&
		Arblib.overlaps(an1.approx, an2.approx)
end

function inv(an::AlgebraicNumber)
	new_minpoly = reverse(an.minpoly)
	new_approx = certify_approx(new_minpoly, Arblib.inv!(Acb(0), an.approx))
	return AlgebraicNumber(new_minpoly, new_approx)
end

# interleave each elemnet of a with n zeros
function deletelastzeros(v)
	i = findlast(!=(0), v)
	return v[begin:i]
end

interleavezeros(v,n) = deletelastzeros(vec(vcat(v', zeros(eltype(v), n, length(v)))))

function root(an::AlgebraicNumber, n::Int)
	out = an
	if n == 0
		# Follow julia's convention for: an^(1/0)=an^(infinity)
		out = AlgebraicNumber([zero(fmpz)], Acb(big(an.approx)^Inf))
	else
		if n < 0
			out = root(inv(an), -n)
		else
			num = Arblib.root(an.approx, n)
			# Need simplify for abs(sqrt(AlgebraicNumber(-1))) == 1
			new_minpoly = simplify(interleavezeros(an.minpoly, n - 1), num)
			new_approx = certify_approx(new_minpoly, num)
			return AlgebraicNumber(new_minpoly, new_approx)
		end
	end
end

root(an, n) = root(AlgebraicNumber(an), Int(n))


sqrt(an::AlgebraicNumber) = root(an, 2)
cbrt(an::AlgebraicNumber) = root(an, 3)

# TODO: special, more efficient cases for ^2 and ^3
function pow2(an::AlgebraicNumber)
	coeffs = an.minpoly
	# first check if it is already in the form of a square root.
	if all(coeffs[2:2:end] .== 0)
		new_minpoly = coeffs[1:2:end]
	else
		new_minpoly = poly_product_from_coeff(
			coeffs,
			(-1) * minus_minpoly(coeffs))[1:2:end]
	end
	new_approx = certify_approx(new_minpoly, an.approx*an.approx)
	return AlgebraicNumber(new_minpoly, new_approx)
end

# multiplication
function *(an1::AlgebraicNumber, an2::AlgebraicNumber)
	out = zero(AlgebraicNumber)
	if !(an1 == 0) && !(an2 == 0)
		if an1 == an2
			out = pow2(an1)
		else
			num = an1.approx * an2.approx
			new_minpoly = simplify(composed_product(an1.minpoly, an2.minpoly), num)
			new_approx = certify_approx(new_minpoly, num)
			out = AlgebraicNumber(new_minpoly, new_approx)
		end
	end
	return out
end

function +(an1::AlgebraicNumber, an2::AlgebraicNumber)
	num = an1.approx + an2.approx
	new_minpoly = simplify(composed_sum(an1.minpoly, an2.minpoly), num)
	new_approx = certify_approx(new_minpoly, num)
	return AlgebraicNumber(new_minpoly, new_approx)
end

# +(an1::AlgebraicNumber{T}, an2::AlgebraicNumber{S}) where {T,S} = +(promote(an1, an2))

function -(an::AlgebraicNumber)
	new_minpoly = minus_minpoly(an.minpoly)
	new_approx = certify_approx(new_minpoly, -an.approx)
	return AlgebraicNumber(new_minpoly, new_approx)
end

-(an1::AlgebraicNumber,an2::AlgebraicNumber) = an1 + (-an2)
/(an1::AlgebraicNumber,an2::AlgebraicNumber) = an1 * (inv(an2))

# the complex conjugate of an algebraic number has the same minimal polynomial
conj(an::AlgebraicNumber) = AlgebraicNumber(an.minpoly, conj(an.approx))
abs(an::AlgebraicNumber) = sqrt(an * conj(an))

real(an::AlgebraicNumber) = (an + conj(an)) * inv(AlgebraicNumber(2))
imag(an::AlgebraicNumber) = (an - conj(an)) * inv(AlgebraicNumber(Complex(0,2)))

# zero(::Type{AlgebraicNumber{T,F}}) where {T,F} = AlgebraicNumber{T,F}(T[0, 1], Complex{F}(0.0), F(Inf))
zero(::Type{AlgebraicNumber}) = AlgebraicNumber([0, 1], 0)
zero(x::AlgebraicNumber) = AlgebraicNumber([0, 1], 0)

one(::Type{AlgebraicNumber}) = AlgebraicNumber([-1, 1], 1)
one(x::AlgebraicNumber) =  AlgebraicNumber([-1, 1], 1)

# TODO
# take roots of a polynomial,
# and return them as algebraic numbers
# function alg_roots(coeff::Vector{Integer})
# end

# compute exp(pi*i*a),
# which is algebraic if a is rational.
function exp_alg(a::IntOrRat)
	# first, obtain polynomial
	minpoly = interleavezeros([-1,1], 2 * denominator(a) - 1)
	# now, select root.
	approx = Arblib.exp_pi_i!(Acb(0),Acb(a))
	# Finally, return minimal polynomial w.r.t. that root
	return AlgebraicNumber(minpoly, approx)
end

exp_alg(an::AlgebraicNumber) = exp_alg(convert(Rational, an))

#using duck catch by exp_alg
cos_alg(a) = real(exp_alg(a))
sin_alg(a) = imag(exp_alg(a))
