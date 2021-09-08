using LinearAlgebra, ForwardDiff, StaticArrays

const J = Diagonal(SA[1,2,3])

function dyno(ω)
    return J\(-cross(ω,(J*ω)))
end

function dynojac(ω)
    return vec(ForwardDiff.jacobian(dyno,ω))
end

function jvp(λ,ω)
    λ'*ForwardDiff.jacobian(dyno,ω)
end
# ForwardDiff.jacobian(dyno,randn(3))
#
xs = randn(3)
λ = randn(3)
# t1 = λ'*ForwardDiff.jacobian(dynojac,xs)
jvp(λ,xs)
