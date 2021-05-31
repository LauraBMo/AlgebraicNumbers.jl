module AlgebraicNumbers

using Nemo
import PolynomialRoots:roots

export AlgebraicNumber
export *,+,-,/,^,root,==,inv
export sqrt,cbrt
export exp_alg,cos_alg,sin_alg
export pow2


include("algebraic.jl")
include("promote.jl")
include("newton.jl")
include("Nemo.jl")
include("Base.jl")

end
