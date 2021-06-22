
@testset "Promote type" begin
    for T in [Int, BigInt]
        @test promote_type(T, AlgebraicNumber) == AlgebraicNumber
        @test promote_type(Rational{T}, AlgebraicNumber) == AlgebraicNumber
        @test promote_type(Complex{T}, AlgebraicNumber) == AlgebraicNumber
        @test promote_type(Complex{Rational{T}}, AlgebraicNumber) == AlgebraicNumber
    end
end

@testset "issue #5" begin
    # testcase of issue #5
    @test AlgebraicNumber(1) + sqrt(AlgebraicNumber(-1)) != AlgebraicNumber(2)
end

@testset "Show" begin
    an = sqrt(AlgebraicNumber(-1)) + 1
    io = IOBuffer()
    show(io, an)
    @test String(take!(io)) == "≃$(Arblib.string_nice(an.approx, 20))"
    # @test String(take!(io)) == "\"≈$(Arblib.string_nice(an.approx, 20))\""
    close(io)
end
