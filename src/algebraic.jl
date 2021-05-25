# Exact representation of algebraic numbers
# (Numbers that are roots of polynomials with integer coefficients)
# And also arithmetic on algebraic numbers,
# including +, -, *, /, and radicals.


import Base.zero,Base.one
import Base.+,Base.-,Base.*,Base./,Base.inv
import Base.abs,Base.conj
import Base.real,Base.imag
import Base.==,Base.hash,Base.show

# see: http://nemocas.org/nemo-0.4.pdf

# An algebraic number,
# consisting of the minimal polynomial of the number,
# an arbitrary-precision approximation of the number,
# and prec which specifies the minimal distance between roots of p
# TODO: apprx has to be complex.
"""
		AlgebraicNumber{T}

An algebraic number consisting of
   - the minimal polynomial of the number: Vector of integers.
   - Complex float point aproximation
   - minimal distance between two roots of minimal polynomial
"""
struct AlgebraicNumber{T <: Integer,F <: AbstractFloat} <: Number
	coeff::Vector{T}
	apprx::Complex{F}
	prec::F
end

# algebraic number from just poly and approximation.
# computes precision and simplifies as well.
function AlgebraicNumber(coeff::Vector{T}, num::S, ::Type{F}=BigFloat) where {T <: Integer,S <: Number,F <: AbstractFloat}
	minpoly = get_minpoly(coeff, num)
	apprx = Complex{F}(num)
	# multiply by 0.5 safety factor
	prec = convert(F, 0.5 * min_pairwise_dist(prec_roots(minpoly)))
	return AlgebraicNumber{T,F}(minpoly, apprx, prec)
end
end

"""
	 AlgebraicNumber(x::T) where {T <: Integer}

Algebraic number from integer.
"""
AlgebraicNumber(x::T) where {T <: Integer} = AlgebraicNumber(BigInt[-x,one(T)], x)

"""
	 AlgebraicNumber(x::T) where {T <: Rational}

Algebraic number from rational.
"""
AlgebraicNumber(x::Rational) = AlgebraicNumber(BigInt[-numerator(x), denominator(x)], x)

"""
	 AlgebraicNumber(x::Complex{T}) where T <: Integer

Algebraic number from ZZ[im].
"""
function AlgebraicNumber(x::Complex{T}) where T <: Integer
	if imag(x) == 0
		return AlgebraicNumber(real(x))
	else
		return AlgebraicNumber(BigInt[imag(x)^2 + real(x)^2,-2 * real(x),one(T)], x)
	end
end

"""
	 AlgebraicNumber(x::Complex{Rational})

Algebraic number from QQ[im].
"""
function AlgebraicNumber(x::Complex{Rational})
	if imag(x) == 0
		return AlgebraicNumber(real(x))
	else
		v = [imag(x)^2 + real(x)^2,-2 * real(x),one(T)]
		return AlgebraicNumber(BigInt.(lcm(denominator.(v)) .* v), x)
	end
end

function poly_from_coeff(a)
	R, x = PolynomialRing(Nemo.FlintZZ, "x")
	sum([a[i] * x^(i - 1) for i = 1:length(a)])
end

Base.hash(an::AlgebraicNumber, h::UInt) = hash(an.coeff, h)

# TODO: only show up to precision
function show(io::IO, an::AlgebraicNumber)
	print(io, "≈")
	# ndigits = max(10, round(Int,ceil(convert(Float64,log(an.prec)/log(10)))))
	show(io, convert(Complex{Float64}, an.apprx))
	# print(io,"...")
end

# get_coeffs(p::Nemo.fmpz_poly) = pointer_to_array(convert(Ptr{Int64}, p.coeffs), (p.length,))
get_coeffs(p::Nemo.fmpz_poly) = [BigInt(Nemo.coeff(p, i)) for i = 0:Nemo.degree(p)]
prec_roots(a::Vector{T}) where {T <: Integer} = PolynomialRoots.roots(convert(Array{BigFloat}, a))
prec_roots(a::PolyElem) = prec_roots(get_coeffs(a))
# TODO: make sure roots returns distinct roots

# Given an algebraic number, find minimum precision required
# to specify it among roots of an.p
# TODO: handle case of repeated roots precisely
function min_pairwise_dist(v::Vector)
	mindist = Inf
	# nontrivial case
	if length(v) > 1
		# find minimum pairwise distance between pairs;
		mindist = minimum(distances(v))
	end
	return mindist
end

# compute distances
distances(x::Number, v::Vector) = abs.(x .- v)

# v of length at least 2
distances(v::Vector, n::Int=length(v)) =
	reduce(vcat, [distances(v[i], v[(i + 1):end]) for i in 1:(n - 1)])

# simplify an algebraic number by reducing p to the minimal polynomial.
function get_minpoly(poly, num)
	# for all factors of an.p, find the one that matches our roots
	minpoly = poly
	# Unless nonlinear polynomial, then already irreducible.
	if length(poly) > 2 # polynomial not linear: factor out
		R, x = PolynomialRing(Nemo.FlintZZ, "x")
		factorisation = collect(Nemo.factor(R(poly)))
		# first, trivial case (one factor)
		if length(factorisation) == 1
			factor = factorisation[1]
		else
			# case where more than one factor exists
			mindists = [minimum(distances(num, roots))
						for roots in prec_roots.(first.(factorisation))]
			(mindist, i) = findmin(mindists)
			factor = factorisation[i]
		end
		minpoly = get_coeffs(first(factor))
	end
	return minpoly
end

moniccoeffs(an::AlgebraicNumber) = an.coeff[begin:(end - 1)] .// an.coeff[end]

function ==(an1::AlgebraicNumber, an2::AlgebraicNumber)
	moniccoeffs(an1) == moniccoeffs(an2) || return false
	return abs(an1.apprx - an2.apprx) < min(an1.prec, an2.prec)
end

inv(an::AlgebraicNumber) = AlgebraicNumber(reverse(an.coeff), inv(an.apprx))

# interleave each elemnet of a with n zeros
interleave(a,n) =  vec(vcat(a', zeros(Int64, n, length(a))))
function root(an::AlgebraicNumber, n::Int64)
	if n == 0
		throw(ArgumentError("n must be nonzero"))
	end
	if n == 1
		return an
	end
	if n < 0
		an = inv(an)
		n = -n
	end
	# TODO: quickly calculate precision
	return AlgebraicNumber(interleave(an.coeff, n - 1), an.apprx^(1 / n))
end

import Base.sqrt
import Base.cbrt
sqrt(an::AlgebraicNumber) = root(an, 2)
cbrt(an::AlgebraicNumber) = root(an, 3)

# TODO: special, more efficient cases for ^2 and ^3
function pow2(an::AlgebraicNumber)
	cfs = an.coeff
	# first check if it is already in the form of a square root.
	if all(cfs[2:2:end] .== 0)
		pp_cfs = cfs
	else
		cfs2 = [iseven(i) ? -cfs[i] : cfs[i] for i = 1:length(cfs)]
		pp = poly_from_coeff(cfs) * poly_from_coeff(cfs2)
		pp_cfs = get_coeffs(pp)
	end
	p2 = pp_cfs[1:2:end]
	return AlgebraicNumber(p2, an.apprx * an.apprx)
end


# partially simplify a polynomial b
# eliminating repeated factors
reduce_repeated_factors(p::Nemo.fmpz_poly) = prod(keys(Nemo.factor(p)))

# multiplication
function *(an1::AlgebraicNumber, an2::AlgebraicNumber)
	if an1 == 0 || an2 == 0
		# TODO: don't handle this explicitly
		return zero(AlgebraicNumber)
	end
	# check if p==q, if then use a more optimized and correct routine
	# if an1.coeff == an2.coeff
	#	return
	# end
	p = composed_product(an1.coeff, an2.coeff)
	return AlgebraicNumber(p, an1.apprx * an2.apprx)
end

function +(an1::AlgebraicNumber, an2::AlgebraicNumber)
	p = composed_sum(an1.coeff, an2.coeff)
	return AlgebraicNumber(p, an1.apprx + an2.apprx)
end

function -(an1::AlgebraicNumber)
	cfs = copy(an1.coeff)
	for i = 1:2:length(cfs)
		cfs[i] = -cfs[i]
	end
	return AlgebraicNumber(cfs, -an1.apprx, an1.prec)
end

-(an1::AlgebraicNumber,an2::AlgebraicNumber) = an1 + (-an2)
/(an1::AlgebraicNumber,an2::AlgebraicNumber) = an1 * (inv(an2))

# the complex conjugate of an algebraic number has the same minimal polynomial
conj(an::AlgebraicNumber) = AlgebraicNumber(an.coeff, conj(an.apprx), an.prec)
abs(an::AlgebraicNumber) = sqrt(an * conj(an))

zero(::Type{AlgebraicNumber}) = AlgebraicNumber(BigInt[0, 1], Complex{BigFloat}(0.0), BigFloat(1.0))
one(::Type{AlgebraicNumber})  = AlgebraicNumber(BigInt[-1,1], Complex{BigFloat}(1.0), BigFloat(1.0))

real(an::AlgebraicNumber) = (an + conj(an)) * AlgebraicNumber(BigInt[1,-2], BigFloat(0.5) + 0im, BigFloat(0.5))
imag(an::AlgebraicNumber) = (an - conj(an)) * AlgebraicNumber(BigInt[1,0,4], BigFloat(-0.5) * im, BigFloat(0.5))

# take roots of a polynomial,
# and return them as algebraic numbers
function alg_roots(coeff::Vector{Integer})
	# TODO
end

confirm_algnumber(b) = sum(b.coeff .* [b.apprx^(i - 1) for i = 1:length(b.coeff)])

# compute exp(pi*i*a),
# which is algebraic if a is rational.
function exp_alg(a::Rational)
	# first, obtain polynomial
	p = interleave(BigInt[-1,1], 2 * denominator(a) - 1)
	# now, select root.
	apprx = exp(im * BigFloat(pi) * a)
	# Finally, return minimal polynomial w.r.t. that root
	return AlgebraicNumber(p, apprx)
end

cos_alg(a::Rational) = real(exp_alg(a))
sin_alg(a::Rational) = imag(exp_alg(a))
