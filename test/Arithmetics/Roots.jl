
@testset "Square root" begin
    for n in randint()
        @test (sqrt(AlgebraicNumber(n)))^2 == n
    end
end

@testset "Cubic root" begin
    for n in randint()
        @test (cbrt(AlgebraicNumber(n)))^3 == n
    end
end

@testset "kth roots" begin
    for n in randint(;length=7), k in 3:15
        an = root(AlgebraicNumber(n), k)^k
        # print("$(an.minpoly)\n")
# print("$(an.approx)\n")
        @test an == n
    end
end
