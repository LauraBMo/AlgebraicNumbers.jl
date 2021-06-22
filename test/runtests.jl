using Test
using AlgebraicNumbers
using Arblib

randint(::Type{T}=BigInt; range=100, length=50, init=-range) where {T}  = T.(rand(init:range, length))
randrat(::Type{T}=BigInt; range=100, length=50) where {T} = randint(T, range=range, length=length) .// randint(T, range=range, length=length, init=1)

isapproxzero = isapprox(0; atol=1e-20)

include("AGtests.jl")

# TODO: add some more tests
