
Base.big(z::Acb) = Complex(BigFloat(Arblib.realref(z)), BigFloat(Arblib.imagref(z)))
Base.big(n::fmpz) = BigInt(n)

certify_approx(poly::Vector{fmpz}, approx::Acb) = certify_approx(big.(poly), big.(approx))

function certify_approx(poly::Vector{T}, approx::S) where {T <: Integer,S <: Number}
    out = Acb(0)
    if !(poly ==  T[0,1]) # problem with certification approx
        t = Variable(:t)
        F = System([vec_to_poly(poly, t)], variables=[t])
        out = certify_approx(F, approx)
    end
    return out
end

function certify_approx(F::System, approx::S) where {S <: Number}
    # Maybe use HomotopyContinuation.certify_solution
    certificate = first(certificates(certify(F, [approx]; show_progress=false)))
    return is_certified(certificate) ?
        first(certified_solution_interval(certificate)) :
        throw("Error certifying solution!")
end
