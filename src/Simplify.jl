
# simplify an algebraic number by reducing p to the minimal polynomial.
function simplify(poly::Vector, num)
    R, x = PolynomialRing(Nemo.FlintZZ, "x")
    factors = Nemo.coeffs.(first.(Nemo.factor(R(poly))))
    val, i = findmin(map(p -> abs(vec_to_poly(BigInt.(p), num)), factors))
    # factors_roots = prec_roots.(factors)
    # # for all factors of an.p, find the one that matches our root
    # (mindist, i) = findmin(minimum.(distances(num).(factors_roots)))
    return collect(factors[i])
end

simplify(poly::Vector, num::Acb) = simplify(poly, big(num))
