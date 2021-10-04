using LinearAlgebra, StaticArrays, BenchmarkTools
using Symbolics

# import Attitude.hat
import Attitude.pdot_from_w
@inline function skew(v::SVector{3, Float64})
    return @SArray [0 -v[3] v[2]; v[3] 0 -v[1]; -v[2] v[1] 0]
end
# @inline function pdot_from_w(p::SVector{3, Float64},ω::SVector{3, Float64})::SVector{3, Float64}
#     return ((1+norm(p)^2)/4) *(   I + 2*(skew(p)^2 + skew(p))/(1+norm(p)^2)   )*ω
# end
function dynamics(model::SAT,x,u)
    p = x[SA[1,2,3]]
    ω = x[SA[4,5,6]]
    pdot = ((1+norm(p)^2)/4) *(   I + 2*(skew(p)^2 + skew(p))/(1+norm(p)^2)   )*ω
    # pdot = pdot_from_w(p,ω)
    α = model.invJ*(u - cross(ω,model.J*ω))
    return SA[pdot[1],pdot[2],pdot[3],α[1],α[2],α[3]]
end

struct SAT
    J::Diagonal{Float64, SVector{3, Float64}}
    invJ::SMatrix{3, 3, Float64, 9}
end


let
    J = Diagonal(SA[1,2,3.0])
    invJ = SMatrix{3,3}(inv(J))

    sat = SAT(J,invJ)

    x = @SArray randn(6)
    u = @SArray randn(3)

    @btime dynamics($sat,$x,$u)
    # @btime hat($u)
    #
    # p = @SArray randn(3)
    # ω = @SArray randn(3)
    #
    # @btime pdot_from_w($p,$ω)
end
