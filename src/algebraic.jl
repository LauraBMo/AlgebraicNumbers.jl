# Exact representation of algebraic numbers
# (Numbers that are roots of polynomials with integer coefficients)
# And also arithmetic on algebraic numbers,
# including +, -, *, /, and radicals.

import Base.zero,Base.one
import Base.+,Base.-,Base.*,Base./,Base.inv
import Base.abs,Base.conj
import Base.real,Base.imag
import Base.==,Base.hash,Base.show
import Base.sqrt,Base.cbrt

# see: http://nemocas.org/nemo-0.4.pdf

floattype(an::AlgebraicNumber{T,F}) where {T,F} = F
complextype(an::AlgebraicNumber{T,F}) where {T,F} = Complex{F}
inttype(an::AlgebraicNumber{T,F}) where {T,F} = T

# algebraic number from just a polynomial and a number approximating one of its roots.
# computes precision and simplifies as well.
function AlgebraicNumber(coeff::Vector{T}, num::S, ::Type{F}=BigFloat) where {T <: Integer,S <: Number,F <: AbstractFloat}
	minpoly, mindist, roots = get_minpoly(coeff, num)
	an = set_an(F, minpoly, num, roots)
	# iswelldefined(mindist, an.prec)
	return an
end

AlgebraicNumber(coeff::Vector{T}, num::Complex{F}) where {T <: Integer,F <: AbstractFloat} = AlgebraicNumber(coeff, num, F)

# TODO check that the new algebraic number is well defined. Use iswelldefined?
# TODO add certificate that num approximates a root of minpoly, see HomotopyContinuation.jl certificate
function set_an(::Type{F}, minpoly::Vector{T}, num::S, roots=prec_roots(minpoly)) where {T <: Integer,S <: Number,F <: AbstractFloat}
	# multiply by 0.3 safety factor (the maximal factor is 1/3, bigger factors do not guarantee ==)
	prec = convert(F, 0.3 * min_pairwise_dist(roots))
	apprx = Complex{F}(num)
	return AlgebraicNumber{T,F}(minpoly, apprx, prec)
end

"""
	 AlgebraicNumber(x::T) where {T <: Integer}

Algebraic number from integer.
"""
AlgebraicNumber(x::T, ::Type{F}=BigFloat) where {T <: Integer,F <: AbstractFloat} =
	set_an(F, [-x, one(x)], x, [x])

"""
	 AlgebraicNumber(x::T) where {T <: Rational}

Algebraic number from rational.
"""
AlgebraicNumber(x::Rational, ::Type{F}=BigFloat) where {F <: AbstractFloat} =
	set_an(F, [-numerator(x), denominator(x)], x, [x])

"""
	 AlgebraicNumber(x::Complex{T}) where T <: Integer

Algebraic number from ZZ[im].
"""
function AlgebraicNumber(x::Complex{T}, ::Type{F}=BigFloat) where {T <: Integer,F <: AbstractFloat}
	if imag(x) == 0
		return AlgebraicNumber(real(x), F)
	else
		return set_an(F, [imag(x)^2 + real(x)^2, -2 * real(x), one(T)], x, [x, conj(x)])
	end
end

"""
	 AlgebraicNumber(x::Complex{Rational})

Algebraic number from QQ[im].
"""
function AlgebraicNumber(x::Complex{Rational{T}}, ::Type{F}=BigFloat) where {T <: Integer,F <: AbstractFloat}
	if imag(x) == 0
		return AlgebraicNumber(real(x), F)
	else
		rat_coeffs = [imag(x)^2 + real(x)^2, -2 * real(x), one(T)]
		return set_an(F, convert_intcoeffs(T, rat_coeffs), x, [x, conj(x)])
	end
end

# get_coeffs(p::Nemo.fmpz_poly) = pointer_to_array(convert(Ptr{Int64}, p.coeffs), (p.length,))
get_coeffs(p::PolyElem, ::Type{T}=BigInt) where {T}  = T.(Nemo.coeffs(p))
prec_roots(a::Vector) = unique(roots(BigFloat.(a)))
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
distances(x::Number) = v -> distances(x, v)
# distance between a polynomial and a number
distance(coeff::Vector{T}, x::Number) where {T <: Integer} = minimum(distances(x, prec_roots(coeff)))

# v of length at least 2
distances(v::Vector) =
	reduce(vcat, [distances(v[i], v[(i + 1):end]) for i in 1:(length(v) - 1)])

# simplify an algebraic number by reducing p to the minimal polynomial.
function get_minpoly(poly::Vector, num)
	if degree(poly) == 1 # polynomial linear
		minpoly = poly
		roots = prec_roots(poly)
		mindist = minimum(distances(num, roots))
	else
		R, x = PolynomialRing(Nemo.FlintZZ, "x")
		factors = first.(Nemo.factor(R(poly)))
		factors_roots = prec_roots.(factors)
		# for all factors of an.p, find the one that matches our root
		(mindist, i) = findmin(minimum.(distances(num).(factors_roots)))
		minpoly = get_coeffs(factors[i], eltype(poly))
		roots = factors_roots[i]
	end
    return minpoly, mindist, roots
end

function ==(an1::AlgebraicNumber, an2::AlgebraicNumber)
	return moniccoeffs(an1) == moniccoeffs(an2) &&
		abs(an1.apprx - an2.apprx) < min(an1.prec, an2.prec)
end

inv(an::AlgebraicNumber) = AlgebraicNumber(reverse(an.coeffs), inv(an.apprx))

# interleave each elemnet of a with n zeros
interleavezeros(a,n) =  vec(vcat(a', zeros(eltype(a), n, length(a))))

function root(an::AlgebraicNumber, n::Int)
	if n == 0
		# Follow julia's convention for: an^(1/0)=an^(infinity)
		return an.apprx^Inf
	end
	if n == 1
		return an
	end
	if n < 0
		return root(inv(a), -n)
	end
	# TODO: quickly calculate precision
	return AlgebraicNumber(interleavezeros(an.coeffs, n - 1), an.apprx^(1 / n))
end

root(n::Int) = (an::AlgebraicNumber) -> root(an, n)

sqrt(an::AlgebraicNumber) = root(an, 2)
cbrt(an::AlgebraicNumber) = root(an, 3)

# TODO: special, more efficient cases for ^2 and ^3
function pow2(an::AlgebraicNumber)
	cfs = an.coeffs
	# first check if it is already in the form of a square root.
	if all(cfs[2:2:end] .== 0)
		pp_cfs = cfs
	else
		cfs2 = (-1) * minus_minpoly(cfs)
		pp_cfs = poly_product_from_coeff(cfs, cfs2)
	end
	p2 = pp_cfs[1:2:end]
	return AlgebraicNumber(p2, an.apprx * an.apprx)
end

# multiplication
function *(an1::AlgebraicNumber, an2::AlgebraicNumber)
	if an1 == 0 || an2 == 0
		# TODO: don't handle this explicitly
		return zero(promote(an1, an2)[1])
	end
	# check if p==q, if then use a more optimized and correct routine
	# if an1.coeff == an2.coeff
	#	return
	# end
	p = composed_product(an1.coeffs, an2.coeffs)
	return AlgebraicNumber(p, an1.apprx * an2.apprx)
end

function +(an1::AlgebraicNumber, an2::AlgebraicNumber)
	p = composed_sum(an1.coeffs, an2.coeffs)
	return AlgebraicNumber(p, an1.apprx + an2.apprx)
end

function -(an::AlgebraicNumber)
	newcoeff = minus_minpoly(an.coeffs)
	newnum = -an.apprx
	return AlgebraicNumber{inttype(an),floattype(an)}(newcoeff, newnum, an.prec)
end

-(an1::AlgebraicNumber,an2::AlgebraicNumber) = an1 + (-an2)
/(an1::AlgebraicNumber,an2::AlgebraicNumber) = an1 * (inv(an2))

# the complex conjugate of an algebraic number has the same minimal polynomial
conj(an::AlgebraicNumber) = AlgebraicNumber(an.coeffs, conj(an.apprx), an.prec)
abs(an::AlgebraicNumber) = sqrt(an * conj(an))

# zero(::Type{AlgebraicNumber{T,F}}) where {T,F} = AlgebraicNumber{T,F}(T[0, 1], Complex{F}(0.0), F(Inf))
zero(::Type{AlgebraicNumber{T,F}}) where {T,F} = AlgebraicNumber(zero(T), F)
zero(::Type{AlgebraicNumber}) = zero(AlgebraicNumber{BigInt,BigFloat})
# zero(x::AlgebraicNumber) = zero(typeof(x))

one(::Type{AlgebraicNumber{T,F}}) where {T,F} = AlgebraicNumber(one(T), F)
one(::Type{AlgebraicNumber}) = one(AlgebraicNumber{BigInt,BigFloat})
# one(x::AlgebraicNumber) = one(typeof(x))

real(an::AlgebraicNumber) = (an + conj(an)) * inv(AlgebraicNumber(inttype(an)(2)))
imag(an::AlgebraicNumber) = (an - conj(an)) * inv(AlgebraicNumber(inttype(an)(2) * Complex{inttype(an)}(0,1)))

# take roots of a polynomial,
# and return them as algebraic numbers
function alg_roots(coeff::Vector{Integer})
	# TODO
end

# compute exp(pi*i*a),
# which is algebraic if a is rational.
function exp_alg(a::Rational{T}, ::Type{F}=BigFloat) where {T <: Integer,F <: AbstractFloat}
	# first, obtain polynomial
	p = interleavezeros(T[-1,1], 2 * denominator(a) - 1)
	# now, select root.
	apprx = exp(im * F(pi) * a)
	# Finally, return minimal polynomial w.r.t. that root
	return AlgebraicNumber(p, apprx, F)
end

exp_alg(an::AlgebraicNumber) = exp_alg(convert(Rational, an), floattype(an))

cos_alg(a::Rational{T}, ::Type{F}=BigFloat) where {T <: Integer,F <: AbstractFloat} = real(exp_alg(a, F))
sin_alg(a::Rational{T}, ::Type{F}=BigFloat) where {T <: Integer,F <: AbstractFloat} = imag(exp_alg(a, F))

cos_alg(an::AlgebraicNumber) = cos_alg(convert(Rational, an), floattype(an))
sin_alg(an::AlgebraicNumber) = sin_alg(convert(Rational, an), floattype(an))
