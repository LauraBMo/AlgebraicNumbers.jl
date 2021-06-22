
using Nemo
p = BigInt[80,0,0,70]
q = BigInt[-80,0,0,70]
AlgebraicNumbers.moniccoeffs(minus2(p))
AlgebraicNumbers.moniccoeffs(AlgebraicNumbers.minus_minpoly(p))

minus2(p)
AlgebraicNumbers.minus_minpoly(p)

p[2:2:length(p)] = -p[2:2:length(p)]

function minus2(v)
    cfs = v
    cfs[2:2:length(cfs)] = -cfs[2:2:length(cfs)]
    return cfs
end

pmq = AlgebraicNumbers.composed_sum([1,0,0,1], [-1,0,0,1])
pmq = AlgebraicNumbers.composed_sum([1,0,1], [1,0,1])


pmq = AlgebraicNumbers.composed_sum([1,0,0,1], [1,0,0,-1])
rr =ComplexF16.(PolynomialRoots.roots(BigFloat.(pmq)))

r1 = PolynomialRoots.roots(BigFloat[1,0,0,1])
r2 = PolynomialRoots.roots(BigFloat[1,0,0,-1])
ComplexF16.(sumpairs(r1,r2))

sumpairs(v,u) = map(x->+(x...), Iterators.product(v,u))


AlgebraicNumbers.composed_product([-1,2], [1,2])

R, x = Nemo.PolynomialRing(Nemo.FlintQQ, "x")
p = BigInt[1,0,0,1]
q = BigInt[-1,0,0,1]
p = BigInt[1,0,1]
q = BigInt[1,0,1]
# compute newton series
n = AlgebraicNumbers.degree(p) * AlgebraicNumbers.degree(q) + 1
a = AlgebraicNumbers.poly_to_newton(p, n, R, x)
b = AlgebraicNumbers.poly_to_newton(q, n, R, x)

# exp series
ee  = R([1 // factorial(BigInt(i)) for i = 0:n])
eei = R([factorial(BigInt(i)) for i = 0:n])

# multiply newton series and invert
m = mullow(R(AlgebraicNumbers.Hadamard(a, ee)), R(AlgebraicNumbers.Hadamard(b, ee)), n + 1)

AlgebraicNumbers.Hadamard(m, eei)
pq = AlgebraicNumbers.newton_to_poly(AlgebraicNumbers.Hadamard(m, eei))
# convert to integer and return
return convert_intcoeffs(T, pq)
n = (length(p) - 1) * (length(q) - 1) + 1
R, x = Nemo.PolynomialRing(Nemo.FlintQQ, "x")
a = to_newton(p, n, R, x)
b = to_newton(q, n, R, x)

# exp series
ee  = R([Nemo.FlintQQ(1 // factorial(BigInt(i))) for i = 0:n])
eei = R([Nemo.FlintQQ(factorial(BigInt(i))) for i = 0:n])

# multiply newton series and invert
m = mullow(hadm(a, ee, R), hadm(b, ee, R), n + 1)
pq = from_newton(to_array(hadm(m, eei, R)))

# convert to integer and return
return map(numerator, pq * lcm(map(denominator, pq)))
Hadamard_exp(p) = [Rational(c) // factorial(i - 1) for (i, c) in enumerate(Nemo.coeffs(p))]
Hadamard_invexp(p) = [Rational(c) * factorial(i - 1) for (i, c) in enumerate(Nemo.coeffs(p))]
AlgebraicNumbers.Hadamard(a,ee)
AlgebraicNumbers.Hadamard(b,ee)

hb = Hadamard_exp(a)
hb = Hadamard_exp(b)

hbxha = mullow(R(Hadamard_exp(a)), R(Hadamard_exp(b)),n+1)
hbha = R(Hadamard_exp(a))*R(Hadamard_exp(b))

newton_to_poly(Hadamard_invexp(hbxha),10)

AlgebraicNumbers.newton_to_poly(Hadamard_invexp(hbha))

from_newton(Hadamard_invexp(hbxha))

Hadamard_invexp(hbxha)
n

AlgebraicNumbers.composed_sum([-1,0,1],[1,0,-1])
AlgebraicNumbers.minus_minpoly([-1,0,1])

AlgebraicNumbers.mi



# x = root(AlgebraicNumber(-54),7)

# x.minpoly
# alg_im = sqrt(AlgebraicNumber(-1))
# real(x)+ imag(x)*alg_im
# real(x)
# imag(x)

# x - conj(x)
# x
# y=-conj(x)
# z=x+y
# num = big(x.approx+y.approx)
# poly = AlgebraicNumbers.composed_sum(x.minpoly, y.minpoly)


# coeffsp = x.minpoly
# coeffsq = y.minpoly
# # compute newton series
# n = AlgebraicNumbers.degree(coeffsp) * AlgebraicNumbers.degree(coeffsq) + 1
# R, var = Nemo.PolynomialRing(Nemo.FlintQQ, "x")
# coeffsA = AlgebraicNumbers.poly_to_Newton(coeffsp, n, R, var)
# coeffsB = AlgebraicNumbers.poly_to_Newton(coeffsq, n, R, var)

# pA = AlgebraicNumbers.Newton_to_poly(BigInt.(Rational.(coeffsA)), n)
# pB = AlgebraicNumbers.Newton_to_poly(BigInt.(Rational.(coeffsB)), n)

# AlgebraicNumbers.simplify(pA, x.approx)
# AlgebraicNumbers.simplify(pB, y.approx)

# # Multiply newton series
# m = mullow(R(AlgebraicNumbers.Hadamard_exp(coeffsA)), R(AlgebraicNumbers.Hadamard_exp(coeffsB)), n + 1)
# # Convert back to polynomial
# pq = AlgebraicNumbers.Newton_to_poly(Rational.(AlgebraicNumbers.Hadamard_invexp(m)), n)

# # convert to integer and return
# vpq = convert(Vector{BigInt}, AlgebraicNumbers.RatVec(pq))

# AlgebraicNumbers.vec_to_poly(vpq, x.approx+y.approx)

# using Nemo
# R, x = PolynomialRing(Nemo.FlintZZ, "x")
# factors = Nemo.coeffs.(first.(Nemo.factor(R(poly))))
# val, i = findmin(map(p -> abs(AlgebraicNumbers.vec_to_poly(BigInt.(p), num)), factors))
# # factors_roots = prec_roots.(factors)
# # # for all factors of an.p, find the one that matches our root
# # (mindist, i) = findmin(minimum.(distances(num).(factors_roots)))
# return integertype(poly).(factors[i])
# new_minpoly = AlgebraicNumbers.simplify(AlgebraicNumbers.composed_sum(x.minpoly, y.minpoly), num)

# new_approx = certify_approx(new_minpoly, num)
# z.minpoly
# x

# 65^(1/7)
# y = Arblib.root(Acb(67),7)

# Arblib.overlaps(x.approx,y)

# BigInt[65, 0, 0, 0, 0, 0, 0, 1]



# using HomotopyContinuation
# num = big(0)
# a = AlgebraicNumber(4)
# b = AlgebraicNumber(-4)
# new_minpoly = AlgebraicNumbers.simplify(AlgebraicNumbers.composed_sum(a.minpoly, b.minpoly), num)
# new_approx = AlgebraicNumbers.certify_approx(new_minpoly, num)

# t = Variable(:t)
# F = System([AlgebraicNumbers.vec_to_poly(new_minpoly, t)], variables=[t])

# # Maybe use HomotopyContinuation.certify_solution

# certify(F, [num]; show_progress=false)
# certificate = first(certificates(certify(F, [num]; show_progress=false)))
# first(certified_solution_interval(certificate))
