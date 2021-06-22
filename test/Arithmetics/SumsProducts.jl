
@testset "Sum" begin
    N = randint()
    M = randint()
    for (n, m) in zip(N, M)
        # print("Suming $n, $m\n")
        @test AlgebraicNumber(n) + AlgebraicNumber(m) == AlgebraicNumber(n + m)
    end
end

@testset "Product" begin
    N = randint()
    M = randint()
    for (n, m) in zip(N, M)
        # print("Multiplying $n, $m\n")
        @test AlgebraicNumber(n) * AlgebraicNumber(m) == AlgebraicNumber(n * m)
    end
end

@testset "Golden ratio" begin
    ϕ = 1 // 2 + sqrt(AlgebraicNumber(5 // 4))
    # As we all know, this has the property that
    @test 1 + 1 / ϕ == ϕ
end

@testset "More products" begin
    for (n, m) in zip(randrat(), randrat())
        @test (sqrt(AlgebraicNumber(n)) * sqrt(AlgebraicNumber(m)))^2 == n * m
        @test (cbrt(AlgebraicNumber(n)) * cbrt(AlgebraicNumber(m)))^3 == n * m
    end
end

@testset "Power of two" begin
    for n in randrat()
        @test AlgebraicNumber(n)^2 == n^2
    @test sqrt(AlgebraicNumbers.pow2(AlgebraicNumber(n))) == abs(n)
    end
end
