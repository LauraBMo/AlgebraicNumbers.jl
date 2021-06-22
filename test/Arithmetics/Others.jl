
@testset "Plastic const" begin
    # see http://mathworld.wolfram.com/PlasticConstant.html
    a = sqrt(AlgebraicNumber(69))
    n = cbrt(9 - a) + cbrt(9 + a)
    p = n * inv(cbrt(AlgebraicNumber(18)))
    @test p - 1 == 1 / (p^4)
    @test p + 1 == p^3
end

@testset "Stackoverflow" begin
    # http://math.stackexchange.com/questions/422233/how-to-find-a-minimal-polynomial-field-theory
    n = sqrt(AlgebraicNumber(9 * 5)) - sqrt(AlgebraicNumber(4 * 7)) + sqrt(AlgebraicNumber(35))
    d = 1 - sqrt(AlgebraicNumber(5)) + sqrt(AlgebraicNumber(7))
    α = n / d
    @test AlgebraicNumbers.getminpoly(α) == [3596, 2312, -280, -156, 19]
    @test Arblib.overlaps(AlgebraicNumbers.getapprox(α), α.approx)
end

@testset "Absolute val" begin
    alg_im = sqrt(AlgebraicNumber(-1))
    @test conj(alg_im) == -alg_im
    @test abs(alg_im) == 1
    for n in randrat()
    @test abs(AlgebraicNumber(n)) == abs(n)
    end
end

@testset "Sin and Cos" begin
    step = 5
    for x in Int.(0:step * pi)
        @test abs(sin_alg(x // step)^2 + cos_alg(x // step)^2) == 1
        @test cos_alg(x // step - 1 // 2) == sin_alg(x // step)
    end
    @test cos_alg(2 // 1 - 1 // 2) == sin_alg(2)
end

@testset "Real and Imag" begin
    alg_im = sqrt(AlgebraicNumber(big(-1)))
    for (n, m) in zip(randrat(), randrat()) # , k in 1:15
        # an = root(AlgebraicNumber(n), k)
        an = AlgebraicNumber(n)
        # print("$(an.minpoly)\n")
        # print("$(an.approx)\n")
        @test real(an) + alg_im * imag(an) == an
    end
end
