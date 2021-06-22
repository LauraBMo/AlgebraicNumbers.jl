
@testset "Integers" begin
    for n in randint()
        # print("$n\n")
        algn = AlgebraicNumber(n)
        @test AlgebraicNumbers.confirm_algnumber(algn) == 0
        @test algn == n
        @test convert(Int, algn) == n
    end
end


@testset "Rationals" begin
    for r in randrat()
        # print("$r\n")
        algr = AlgebraicNumber(r)
        @test isapproxzero(AlgebraicNumbers.confirm_algnumber(algr))
        @test algr == r
        @test convert(Rational, algr) == r
    end
end

@testset "Complex Int" begin
    for c in randint() + randint() * Complex(0, 1)
        # print("$c\n")
        algc = AlgebraicNumber(c)
        @test AlgebraicNumbers.confirm_algnumber(algc) == 0
        @test algc == c
        @test convert(Complex, algc) == c
    end
end

@testset "Complex Rat" begin
    for c in randrat() + randrat() * Complex(0, 1)
        # print("$c\n")
        algc = AlgebraicNumber(c)
        @test isapproxzero(AlgebraicNumbers.confirm_algnumber(algc))
        @test algc == c
        @test convert(Complex, algc) == c
    end
end
