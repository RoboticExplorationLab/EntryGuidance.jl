using Convex, Mosek, MosekTools
# using JuMP

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
    # @infiltrate
    # error()
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
    # γ = 0.0
    α = 1e-8/length(U)
    # α = 0.0
    # β = 1e-8/length(U)
    σ_0 = [X[i][7] for i = 1:length(X)]
    #                        miss distance                                   bank angle                          regularizer on steps
    # @infiltrate
    # error()
    # problem = minimize( 1000*sumsquares( ( Qn*(X[N][1:3] + δx[1:3,N]) - xf[1:3])  )+β*sumsquares(vec(mat_from_vec(U)) + δu), cons)
    problem = minimize( 1000*sumsquares( ( Qn*(X[N][1:3] + δx[1:3,N]) - xf[1:3])  ) + γ*sumsquares( σ_0 + vec(δx[7,:]) ) + α*sumsquares(vec(δu)), cons)
    # problem = minimize( 1000*sumsquares( ( Qn*(X[N][1:3] + δx[1:3,N]) - xf[1:3])  ) + γ*norm( σ_0 + vec(δx[7,:]) ) + α*sumsquares(vec(δu)), cons)
    Convex.solve!(problem, () -> Mosek.Optimizer())

    cX = vec_from_mat(mat_from_vec(X) + evaluate(δx))
    # @infiltrate
    # error()
    # cU = vec_from_mat(mat_from_vec(U) + evaluate(δx))
    cU = deepcopy(U)
    # @infiltrate
    # error()
    for i = 1:length(U)
        cU[i][1] += δu.value[i][1]
    end

    return cX, cU, norm(evaluate(δu))
end
# function eg_mpc(model::EntryVehicle,A,B,X,U,xf, mpc_iteration)
#
#     # quick size cehcks
#     @assert (length(A) == length(B))
#     @assert (length(X) == length(U)+1)
#     @assert (length(A) == length(U))
#     N = length(X)
#     x0 = copy(X[1])
#
#     # starting δx
#     δx0 = x0 - X[1]
#     @assert norm(x0 - X[1]) ≈ 0
#
#     # variables
#     δx = Variable(6,N)
#     δu = Variable(2,N-1)
#
#     # dynamics constraints
#     cons = Constraint[ δx[:,1]==δx0 ]
#
#     for i = 1:N-1
#         push!(cons, δx[:,i+1] == sparse(A[i])*δx[:,i] + B[i]*δu[:,i])
#     end
#
#     # lift vector constraints
#     for i = 1:N-1
#         push!(cons, norm(U[i] + δu[:,i]) <= 1.0)
#     end
#
#
#
#     rr = normalize(xf[1:3])
#     Qn = I - rr*rr'
#
#     # trust region
#     # push!(cons, norm( Qn*δx[1:3,N]) <= 5)
#     # push!(cons,norm( Qn*( (X[N][1:3] + δx[1:3,N]) - xf[1:3])) <= norm( Qn*( (X[N][1:3]) - xf[1:3])  ) )
#
#     γ = 10000         # miss distance penalty
#     α = 1e-4/length(U) # regularizer
#     β = 1/length(U) # control penalty
#
#
#     problem = minimize(γ*sumsquares( Qn*( (X[N][1:3] + δx[1:3,N]) - xf[1:3])  ) +  α*sumsquares(vec(δu)) + β*sumsquares(vec(mat_from_vec(U)) + vec(δu)), cons)
#
#
#     # problem = minimize( norm( Qn*( (X[N][1:3] + δx[1:3,N]) - xf[1:3])  ) , cons)
#     Convex.solve!(problem, () -> Mosek.Optimizer())
#
#
#
#     cX = vec_from_mat(mat_from_vec(X) + evaluate(δx))
#     cU = vec_from_mat(mat_from_vec(U) + evaluate(δu))
#
#     return cX, cU, norm(evaluate(δu))
#
# end
