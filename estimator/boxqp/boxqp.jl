using LinearAlgebra, OSQP, StaticArrays, SparseArrays
using Infiltrator

function genprob()

    # n = 100;
    # A = randn(n,n); A = A'*A + I;
    # A = Symmetric(A)
    # A2 = copy(A)
    # A3 = copy(A)
    #
    # @btime LAPACK.potrf!('U',$A)

    # @btime cholesky($A2)
    #
    # @btime cholesky($A3)
    n = 5
    H = randn(n,n);H = H'*H
    q = randn(n)

    blo = -ones(n)
    bhi = ones(n)

    m = OSQP.Model()

    OSQP.setup!(m; P=sparse(H), q=q, A=sparse(Matrix(I(n))), l=blo, u=bhi)

    results = OSQP.solve!(m)

    @show results.x

    H = SMatrix{n,n}(H)
    q = SVector{n}(q)
    blo = SVector{n}(blo)
    bhi = SVector{n}(bhi)

    @btime boxqp($H,$q,$blo,$bhi)
    # boxqp(H,q,blo,bhi)

    return nothing
end

function boxqp(H,q,blo,bhi)

    # @btime begin
    n = length(q)

    x = @SVector zeros(5)
    x = clamp.(x,blo,bhi)
    cc = NaN*zeros(5)
    ff = NaN*zeros(5)
    # end
    for i = 1:1
        g = q + H*x
    #
        # find clamped and fre indicees
        c = SA{Int64}[]
        f = SA{Int64}[]
        for ii = 1:n
            # if ((x[ii] == blo[ii]) & (g[ii]>0)) || ((x[ii] == bhi[ii]) & (g[ii]<0))
            if 4>3
                c = push(c,Int(ii))
            else
                f = push(f,Int(ii))
            end
        end
    #
    #     # Δxf = -H[f,f]\g[f]
    #     # Δx = @SVector [k∈f ? Δxf[]]
    #     #
    #     # α = 1.0
    #     # J = 0.5*x'*H*x + q'*x
    #     #
    #     # for j = 1:10
    #     #     # x2 = copy(x)
    #     #     # x2[f] += α*Δxf # this is the current problem
    #     #     x2 = @SVector [ k ∈ f ? x[k] + α*Δx[k] : x[k] for k = 1:5]
    #     #     x2 = clamp.(x2,blo,bhi)
    #     #     if (0.5*x2'*H*x2 + q'*x2 )< J
    #     #         x = copy(x2)
    #     #         break
    #     #     else
    #     #         α *= 0.5
    #     #     end
    #     # end
    #     # if norm(g[f])<1e-5
    #     #     @info "success"
    #     #     break
    #     # end
    #
    #
    #
    end




end

genprob()
