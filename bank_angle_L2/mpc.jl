using Convex, Mosek, MosekTools



function eg_mpc(model::EntryVehicle,A,B,X,U,xf)

    @assert (length(A) == length(B))
    @assert (length(X) == length(U)+1)
    @assert (length(A) == length(U))
    N = length(X)
    x0 = copy(X[1])

    # starting δx
    δx0 = x0 - X[1]
    @assert norm(x0 - X[1]) ≈ 0

    # variables
    δx = Variable(7,N)
    δu = Variable(N-1)

    # dynamics constraints
    cons = Constraint[ δx[:,1]==δx0 ]

    for i = 1:N-1
        push!(cons, δx[:,i+1] == sparse(A[i])*δx[:,i] + B[i]*δu[i])
    end

    # bank angle constraint
    for i = 1:N
        push!(cons, X[i][7] + δx[7,i] <=   3.14159)
        push!(cons, X[i][7] + δx[7,i] >=  -3.14159)
    end

    α = 1e-1
    β = 0e1
    rr = normalize(xf[1:3])
    Qn = I - rr*rr'

    # trust region stuff
    # push!(cons, norm(δx[:,N])<=200.0)
    for i = 1:N
        push!(cons, δx[7,i] <=  deg2rad(20))
        push!(cons, δx[7,i] >= -deg2rad(20))
    end

    # cost function
    γ = 1e-8/length(U)
    α = 0e-3/length(U)
    σ_0 = [X[i][7] for i = 1:length(X)]
    #                        miss distance                                   bank angle                          regularizer on steps
    problem = minimize( norm( ( Qn*(X[N][1:3] + δx[1:3,N]) - xf[1:3])  ) + γ*sumsquares( σ_0 + vec(δx[7,:]) ) + α*sumsquares(vec(δu)), cons)
    Convex.solve!(problem, () -> Mosek.Optimizer())

    cX = vec_from_mat(mat_from_vec(X) + evaluate(δx))
    cU = U + evaluate(δu)

    # bank = [cX[i][7] for i = 1:length(cX)]
    # mat"
    # figure
    # hold on
    # plot($bank)
    # hold off
    # "
    # @infiltrate
    # error()
    return cX, cU
end
