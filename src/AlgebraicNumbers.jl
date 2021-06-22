module AlgebraicNumbers

using Nemo
using Arblib
using HomotopyContinuation
using LinearAlgebra:dot
# using PolynomialRoots:roots

export AlgebraicNumber
export *,+,-,/,^,root,==,inv
export sqrt,cbrt
export exp_alg,cos_alg,sin_alg

include("VectorRational.jl")
include("AlgebraicNumber.jl")
include("Nemo.jl")
include("Base.jl")
include("Newton.jl")
include("ComposedOperations.jl")
include("Simplify.jl")
include("Certify.jl")
include("Arithmetics.jl")

end
