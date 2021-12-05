using QDLDL, SparseArrays, SuiteSparse
using LinearAlgebra
using Infiltrator



# function iterative_ref(F,A,x,b)
#     r = copy(b)
#     x = zeros(length(b))
#     for i = 1:10
#         x += F\r
#         r = b - A*x
#         # @show norm(r)
#         if norm(r)<1e-8
#             break
#         end
#         if i == 10
#             @infiltrate
#             error()
#         end
#     end
#     return x
# end
function iterative_ref(F,A,x,b,r,Δx)
    @. r = b
    for i = 1:10
        @. Δx = r
        solve!(F,Δx)
        @. x += Δx
        mul!(Δx,A,x)
        @. r = b - Δx
        if dot(r,r)<1e-20
            break
        end
    end
    return nothing
end

function ttttt()

    n = 150
    A = sprand(n,n,0.1)

    A = A'*A + I

    F = qdldl(A + 1e-4*I)
    b = randn(n)
    x = zeros(n)
    r = zeros(n)
    Δx = zeros(n)
    @btime x = iterative_ref($F,$A,$x,$b,$r,$Δx)

    @show norm(x - A\b)
    return nothing
end

ttttt()
