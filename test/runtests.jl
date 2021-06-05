using Test
using AlgebraicNumbers

randint(::Type{T}=BigInt; range=100, length=50, init=-range) where {T}  = T.(rand(init:range, length))
randrat(::Type{T}=BigInt; range=100, length=50) where {T} = randint(T, range=range, length=length) .// randint(T, range=range, length=length, init=1)

@testset "Integers" begin
    N = randint()
    algN = AlgebraicNumber.(N)
    @test all(AlgebraicNumbers.confirm_algnumber.(algN) .== 0)
    @test all(algN .== N)
    intN = [convert(Int, n) for n in N]
    @test all(intN .== N)
end

@testset "Rationals" begin
    N = randrat()
    algN = AlgebraicNumber.(N)
    @test all(abs.(AlgebraicNumbers.confirm_algnumber.(algN)) .< 1e-20)
    @test all(algN .== N)
    intN = [convert(Rational, n) for n in N]
    @test all(intN .== N)
end
    
@testset "Complex Int" begin
    N = randint() .+ randint() * Complex{BigInt}(0, 1)
    algN = AlgebraicNumber.(N)
    @test all(AlgebraicNumbers.confirm_algnumber.(algN) .== 0)
    @test all(algN .== N)
    intN = [convert(Complex, n) for n in N]
    @test all(intN .== N)
end

@testset "Complex Rational" begin
    N = randrat() .+ randrat() * Complex{Rational{BigInt}}(0, 1)
    algN = AlgebraicNumber.(N)
    @test all(abs.(AlgebraicNumbers.confirm_algnumber.(algN)) .< 1e-20)
    @test all(algN .== N)
    intN = [convert(Complex, n) for n in N]
    @test all(intN .== N)
end

@testset "Square root" begin
    N = randint()
    sqn = sqrt.(AlgebraicNumber.(N))
    @test all(sqn.^2 .== N)
    intN = [convert(BigInt, r^2) for r in sqn]
    @test all(intN .== N)
end

@testset "Cubic root" begin
    N = randint()
    cbn = cbrt.(AlgebraicNumber.(N))
    @test all(cbn.^3 .== N)
    intN = [convert(BigInt, r^3) for r in cbn]
    @test all(intN .== N)
end

@testset "kth roots" begin
    N = randint(length=10)
    K = randint(Int; range=30, length=10, init=2)
    roots = [root(AlgebraicNumber(n), k) for (n, k) in zip(N, K)]
    @test all([r^k == n for (r, k, n) in zip(roots, K, N)])
    intN = [convert(BigInt, r^k) for (r, k) in zip(roots, K)]
    @test all(intN .== N)
end

@testset "issue #5" begin
    # testcase of issue #5
    @test AlgebraicNumber(1) + sqrt(AlgebraicNumber(-1)) != AlgebraicNumber(2)
end

print("Sums and products:\n")

@testset "Sum" begin
    N = randint()
    M = randint()
    @test all(AlgebraicNumber.(N) .+ AlgebraicNumber.(M) .== AlgebraicNumber.(N .+ M))
    @test all(AlgebraicNumber.(N) .+ AlgebraicNumber.(M) .== AlgebraicNumber.(N .+ M))
end

@testset "Product" begin
    N = randint()
    M = randint()
    @test all(AlgebraicNumber.(N) .* AlgebraicNumber.(M) .== AlgebraicNumber.(N .* M))
    @test all(AlgebraicNumber.(N) .* AlgebraicNumber.(M) .== AlgebraicNumber.(N .* M))
end

    @testset "Golden ratio" begin
    ϕ = 1 // 2 + sqrt(AlgebraicNumber(5 // 4))
    # As we all know, this has the property that
    @test 1 + 1 / ϕ == ϕ
end

@testset "Plastic constant" begin
    # see http://mathworld.wolfram.com/PlasticConstant.html
    a = sqrt(AlgebraicNumber(69))
    n = cbrt(9 - a) + cbrt(9 + a)
    p = n * inv(cbrt(AlgebraicNumber(18)))
    @test p - 1 == 1 / (p^4)
    @test p + 1 == p^3
end

@testset "Stackoverflow example" begin
    # http://math.stackexchange.com/questions/422233/how-to-find-a-minimal-polynomial-field-theory
    n = sqrt(AlgebraicNumber(9 * 5)) - sqrt(AlgebraicNumber(4 * 7)) + sqrt(AlgebraicNumber(35))
d = 1 - sqrt(AlgebraicNumber(5)) + sqrt(AlgebraicNumber(7))
    α = n / d
    @test α.coeffs == BigInt[3596, 2312, -280, -156, 19]
end

@testset "Absolute value" begin
    alg_im = sqrt(AlgebraicNumber(-1))
    @test conj(alg_im) == -alg_im
    @test abs(alg_im) == 1
    N = randrat()
    @test all(abs.(AlgebraicNumber.(N)) == abs.(N))
end

@testset "Real and Imaginary parts" begin
    alg_im = sqrt(AlgebraicNumber(-1))
    N = randint(length=2)
    K = randint(Int; range=5, length=2, init=2)
    roots = [root(AlgebraicNumber(n), k) for (n, k) in zip(N, K)]
    @test all(real.(roots) .+ alg_im * imag.(roots) == roots)
    all(real.(roots) .+ alg_im * imag.(roots) .== roots)
    real.(roots) .+ alg_im * imag.(roots) .== roots
    rts2 = real.(roots) .+ alg_im * imag.(roots)
    x = roots[2]
    y = rts2[2]
    x.coeffs, x.apprx, x.prec
    y.coeffs, y.apprx, y.prec
    AN(an) = an.coeffs, an.apprx, an.prec
    AN(x - conj(x))
    AN(x)
imag(an::AlgebraicNumber) = (an - conj(an)) * inv(AlgebraicNumber(inttype(an)(2) * Complex{inttype(a)}(0, 1)))

r = PolynomialRoots.roots(BigFloat[3596, 2312, -280, -156, 19])
   p1 = AlgebraicNumber(BigInt[3596, 2312, -280, -156, 19], r[1])
   p2 = AlgebraicNumber(BigInt[3596, 2312, -280, -156, 19], r[2])
   p3 = AlgebraicNumber(BigInt[3596, 2312, -280, -156, 19], r[3])
   xx = p1 + p2
    AN(xx)
PolynomialRoots.roots(BigFloat.(xx.coeffs))

    N = randrat() .+ randrat() * Complex{Rational{BigInt}}(0, 1)
    N = abs.(randrat())
    N = cbrt.(AlgebraicNumber.(N))
AN.(N .- conj.(N))
unique(AlgebraicNumbers.degree.(N .+ (-conj.(N))))
AN.(-conj.(N))
    AlgebraicNumbers.composed_sum([80,0,0,70], [-80,0,0,70])
    AN.(N)
PolynomialRoots.roots(BigFloat.(xx.coeffs))
    end

@testset "Sin and Cos" begin
    step = 15
    @test all([abs(sin_alg(x // step)^2 + cos_alg(x // step)^2) == 1 for x in BigInt.(0:step * pi)])
    @test all([cos_alg(x // step - 1 // 2) == sin_alg(x // step) for x in BigInt.(0:step * pi)])
    @test cos_alg(2 // 1 - 1 // 2) == sin_alg(x)
end

@testset "Power of two" begin
    N = abs.(randrat())
    @test all(pow2.(N) == N.^2)
    @test all(sqrt.(pow2.(N)) .== N)
end

@testset "Show" begin
    a = IOBuffer()
    show(a, sqrt(AlgebraicNumber(-1)) + 1)
    @test convert(UTF8String, takebuf_array(a)) == "≈1.0 + 1.0im"
end


# TODO: add some more tests

@testset "More products" begin
    N = randint(init=1)
    M = randint(init=1)
    @test all(sqrt.(AlgebraicNumber.(N)) .* sqrt.(AlgebraicNumber.(M)) == sqrt.(AlgebraicNumber.(N .* M)))
end
