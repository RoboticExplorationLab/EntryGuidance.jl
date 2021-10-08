using Convex, Mosek, MosekTools
# using JuMP

# function convec(model::EntryVehicle,X,U)
#     nx = 7
#     nu = 1
#     N = length(X)
#     c_dyn = zeros()
function eg_mpc2(model::EntryVehicle,A,B,X,U,xf)

    # @show "in"
    # nx = size(B[1],1)
    # nu = size(B[1],2)
    nx = 7
    nu = 1
    N = length(X)

    idx_x = [(i-1)*(nx+nu) .+ (1:nx) for i = 1:length(X)]
    idx_u = [(i-1)*(nx+nu) .+ nx .+ (1:nu) for i = 1:(length(X)-1)]

    nz = (N*nx) + (N-1)*nu
    nc = (N-1)*nx + nx + (N)*1
    A_eq = spzeros(nc,nz)
    # b_eq = zeros(nc)
    # @infiltrate
    # error()
    idx_c = [(i-1)*(nx) .+ (1:nx) for i = 1:(N-1)]
    for i = 1:(N-1)
        @show i
        A_eq[idx_c[i],idx_x[i]]   = A[i]
        A_eq[idx_c[i],idx_u[i]]   = B[i]
        A_eq[idx_c[i],idx_x[i+1]] = -I(nx)
    end
    A_eq[(N-1)*nx .+ (1:nx), idx_x[1]] = I(nx)
    # b_eq[(N-1)*nx .+ (1:nx)] = ones(nx)

    for i = 1:N
        A_eq[N*nx + i,idx_x[i][7]] = 1
    end

    low_eq = zeros((N-1)*nx)
    up_eq = copy(low_eq)
    low_x0 = zeros(nx)
    up_x0 = zeros(nx)
    up_tr = deg2rad(20)*ones(N)
    low_tr = -deg2rad(20)*ones(N)

    A = A_eq
    Lo = [low_eq;low_x0;low_tr]
    Up = [up_eq;up_x0;up_tr]

    R = 1
    P = spzeros(nz,nz)
    q = zeros(nz)
    for i = 1:(N-1)
        P[idx_u[i],idx_u[i]] = [R]
        q[idx_u[i]] = [R*U[i][1]]
    end
    rr = normalize(xf[1:3])
    Qn = 1000*(I - rr*rr')
    P[idx_x[N][1:3],idx_x[N][1:3]] = Qn'*Qn
    q[idx_x[N][1:3]] = -(Qn'*Qn)'*(xf[1:3] - X[N][1:3])

    # solve OSQP
    osqp = OSQP.Model()
    OSQP.setup!(osqp; P=P, q=q, A=A, l=(Lo), u=(Up), eps_abs = 1e-7, eps_rel = 1e-7, polish = 1)
    results = OSQP.solve!(osqp)

    δx = [results.x[idx_x[i]] for i = 1:(N)]
    δu = [results.x[idx_u[i]] for i = 1:(N-1)]

    cX = X + δx
    cU = U + δu

    return cX, cU, norm(results.x)
end
    # @infiltrate
    # error()
    # mat"
    # spy([$A_eq $b_eq])
    # "
    # @infiltrate
    # error()
    #
    #
    # @assert (length(A) == length(B))
    # @assert (length(X) == length(U)+1)
    # @assert (length(A) == length(U))
    # N = length(X)
    # x0 = copy(X[1])
    #
    # # starting δx
    # δx0 = x0 - X[1]
    # @assert norm(x0 - X[1]) ≈ 0
    #
    # # variables
    # δx = Variable(7,N)
    # δu = Variable(N-1)
    #
    # # dynamics constraints
    # # @infiltrate
    # # error()
    # cons = Constraint[ δx[:,1]==δx0 ]
    #
    # for i = 1:N-1
    #     push!(cons, δx[:,i+1] == sparse(A[i])*δx[:,i] + B[i]*δu[i])
    # end
    #
    # # bank angle constraint
    # for i = 1:N
    #     push!(cons, X[i][7] + δx[7,i] <=   3.14159)
    #     push!(cons, X[i][7] + δx[7,i] >=  -3.14159)
    # end
    #
    # α = 1e-1
    # β = 0e1
    # rr = normalize(xf[1:3])
    # Qn = I - rr*rr'
    #
    # # trust region stuff
    # # push!(cons, norm(δx[:,N])<=200.0)
    # for i = 1:N
    #     push!(cons, δx[7,i] <=  deg2rad(20))
    #     push!(cons, δx[7,i] >= -deg2rad(20))
    # end
    #
    # # cost function
    # γ = 1e-8/length(U)
    # # γ = 0.0
    # α = 0e-8/length(U)
    # # α = 0.0
    # β = 1e-16/length(U)
    # σ_0 = [X[i][7] for i = 1:length(X)]
    # #                        miss distance                                   bank angle                          regularizer on steps
    # # @infiltrate
    # # error()
    # # problem = minimize( 1000*sumsquares( ( Qn*(X[N][1:3] + δx[1:3,N]) - xf[1:3])  )+β*sumsquares(vec(mat_from_vec(U)) + δu), cons)
    # problem = minimize( 1000*sumsquares( ( Qn*(X[N][1:3] + δx[1:3,N]) - xf[1:3])  ) + γ*sumsquares( σ_0 + vec(δx[7,:]) ) + α*sumsquares(vec(δu)) + β*sumsquares(vec(mat_from_vec(U)) + δu), cons)
    # # problem = minimize( 1000*sumsquares( ( Qn*(X[N][1:3] + δx[1:3,N]) - xf[1:3])  ) + γ*norm( σ_0 + vec(δx[7,:]) ) + α*sumsquares(vec(δu)), cons)
    # Convex.solve!(problem, () -> Mosek.Optimizer())
    #
    # cX = vec_from_mat(mat_from_vec(X) + evaluate(δx))
    # # @infiltrate
    # # error()
    # # cU = vec_from_mat(mat_from_vec(U) + evaluate(δx))
    # cU = deepcopy(U)
    # # @infiltrate
    # # error()
    # for i = 1:length(U)
    #     cU[i][1] += δu.value[i][1]
    # end

    # return cX, cU, norm(evaluate(δu))
# end
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
    α = 0e-8/length(U)
    # α = 0.0
    β = 1e-16/length(U)
    σ_0 = [X[i][7] for i = 1:length(X)]
    #                        miss distance                                   bank angle                          regularizer on steps
    # @infiltrate
    # error()
    # problem = minimize( 1000*sumsquares( ( Qn*(X[N][1:3] + δx[1:3,N]) - xf[1:3])  )+β*sumsquares(vec(mat_from_vec(U)) + δu), cons)
    problem = minimize( 1000*sumsquares( ( Qn*(X[N][1:3] + δx[1:3,N]) - xf[1:3])  ) + γ*sumsquares( σ_0 + vec(δx[7,:]) ) + α*sumsquares(vec(δu)) + β*sumsquares(vec(mat_from_vec(U)) + δu), cons)
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
