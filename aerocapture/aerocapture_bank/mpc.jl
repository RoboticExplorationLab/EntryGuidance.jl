using Convex, Mosek, MosekTools
# using JuMP
using OSQP


function eg_mpc(model::EntryVehicle,A,B,X,U,ϵ_f)

    # quick size cehcks
    @assert (length(A) == length(B))
    @assert (length(X) == length(U)+1)
    @assert (length(A) == length(U))
    N = length(X)
    x0 = copy(X[1])

    # starting δx
    δx0 = x0 - X[1]
    @assert norm(x0 - X[1]) ≈ 0

    # variables
    δx = Variable(8,N)
    δu = Variable(1,N-1)

    # dynamics constraints
    cons = Constraint[ δx[:,1]==δx0 ]

    for i = 1:N-1
        push!(cons, δx[:,i+1] == sparse(A[i])*δx[:,i] + B[i]*δu[:,i])
    end

    γ = 1e5         # miss distance penalty
    α = 1 # regularizer
    β = 1 # control penalty

    for i = 1:length(X)
        push!(cons, abs(X[i][7] + δx[7,i]) <= deg2rad(180) )
    end

    p = square(X[N][8] + δx[8,N] - ϵ_f)

    # trust region
    push!(cons, norm(δx[7,:],Inf) <= deg2rad(10))
    # for i = 1:N
    #     push!(cons, abs(δx[7,:]) <= deg2rad(10))
    # end

    problem = minimize(γ*p +  α*sumsquares(vec(δu)), cons)

    Convex.solve!(problem, () -> Mosek.Optimizer())

    cX = vec_from_mat(mat_from_vec(X) + evaluate(δx))
    cU = vec_from_mat(mat_from_vec(U) + evaluate(δu))

    # @infiltrate
    return cU

end


function eg_mpc_quad(model::EntryVehicle,A,B,X,U,ϵ_f)

    # sizes for state and control
    nx = 8
    nu = 1
    N = length(X)

    # indicees for state and control
    idx_x = [(i-1)*(nx+nu) .+ (1:nx) for i = 1:length(X)]
    idx_u = [(i-1)*(nx+nu) .+ nx .+ (1:nu) for i = 1:(length(X)-1)]

    # constraint jacobian (in this form L≦Az≤U)
    nz = (N*nx) + (N-1)*nu
    nc = (N-1)*nx + nx #+ (N)*1
    A_eq = spzeros(nc,nz)
    A_ineq = spzeros(N,nz) # state constraints

    # dynamics constraint (equality)
    idx_c = [(i-1)*(nx) .+ (1:nx) for i = 1:(N-1)]
    for i = 1:(N-1)
        A_eq[idx_c[i],idx_x[i]]   = A[i]
        A_eq[idx_c[i],idx_u[i]]   = B[i]
        A_eq[idx_c[i],idx_x[i+1]] = -I(nx)
    end
    A_eq[(N-1)*nx .+ (1:nx), idx_x[1]] = I(nx)

    # state constraints on δσ (inequality)
    for i = 1:N
        A_ineq[i,idx_x[i][7]] = 1
    end

    # constraint bounds
    low_dyneq = zeros((N-1)*nx)
    up_dyneq = copy(low_dyneq)
    low_x0 = zeros(nx)
    up_x0 = zeros(nx)
    up_tr = deg2rad(10)*ones(N)
    low_tr = -deg2rad(10)*ones(N)

    low_eq = [low_dyneq;low_x0]
    up_eq = [up_dyneq;up_x0]
    # stack everything up

    A = [A_eq;A_ineq]
    Lo = [low_eq;low_tr]
    Up = [up_eq;up_tr]

    # cost function terms
    R = 1
    P = spzeros(nz,nz)
    q = zeros(nz)
    for i = 1:(N-1)
        P[idx_u[i],idx_u[i]] = [R]
        q[idx_u[i]] = [R*U[i][1]]
    end
    γ = 1e5
    P[idx_x[N][8],idx_x[N][8]] = γ
    q[idx_x[N][8]] = γ*(X[N][8] - ϵ_f)
    # rr = normalize(xf[1:3])
    # Qn = 1000*(I - rr*rr')
    # P[idx_x[N][1:3],idx_x[N][1:3]] = Qn'*Qn
    # q[idx_x[N][1:3]] = -(Qn'*Qn)'*(xf[1:3] - X[N][1:3])

    # solve using OSQP
    osqp = OSQP.Model()
    OSQP.setup!(osqp; P=P, q=q, A=A, l=(Lo), u=(Up), eps_abs = 1e-8, eps_rel = 1e-8, max_iter = 10000,polish = 1)
    results = OSQP.solve!(osqp)

    # recover solution
    δx = [results.x[idx_x[i]] for i = 1:(N)]
    δu = [results.x[idx_u[i]] for i = 1:(N-1)]

    # apply the δ's
    cX = X + δx
    cU = U + δu

    return cU
end


function eg_mpc_pdip(model::EntryVehicle,A,B,X,U,ϵ_f)

    # sizes for state and control
    nx = 8
    nu = 1
    N = length(X)

    # indicees for state and control
    idx_x = [(i-1)*(nx+nu) .+ (1:nx) for i = 1:length(X)]
    idx_u = [(i-1)*(nx+nu) .+ nx .+ (1:nu) for i = 1:(length(X)-1)]

    # constraint jacobian (in this form L≦Az≤U)
    nz = (N*nx) + (N-1)*nu
    nc = (N-1)*nx + nx #+ (N)*1
    A_eq = spzeros(nc,nz)
    A_ineq = spzeros(N,nz) # state constraints

    # dynamics constraint (equality)
    idx_c = [(i-1)*(nx) .+ (1:nx) for i = 1:(N-1)]
    for i = 1:(N-1)
        A_eq[idx_c[i],idx_x[i]]   = A[i]
        A_eq[idx_c[i],idx_u[i]]   = B[i]
        A_eq[idx_c[i],idx_x[i+1]] = -I(nx)
    end
    A_eq[(N-1)*nx .+ (1:nx), idx_x[1]] = I(nx)

    # state constraints on δσ (inequality)
    for i = 1:N
        A_ineq[i,idx_x[i][7]] = 1
    end

    # constraint bounds
    low_dyneq = zeros((N-1)*nx)
    up_dyneq = copy(low_dyneq)
    low_x0 = zeros(nx)
    up_x0 = zeros(nx)
    up_tr = deg2rad(10)*ones(N)
    low_tr = -deg2rad(10)*ones(N)

    low_eq = [low_dyneq;low_x0]
    up_eq = [up_dyneq;up_x0]
    # stack everything up

    A = [A_eq;A_ineq]
    Lo = [low_eq;low_tr]
    Up = [up_eq;up_tr]

    # cost function terms
    R = 1
    P = spzeros(nz,nz)
    q = zeros(nz)
    for i = 1:(N-1)
        P[idx_u[i],idx_u[i]] = [R]
        q[idx_u[i]] = [R*U[i][1]]
    end
    γ = 1e5
    P[idx_x[N][8],idx_x[N][8]] = γ
    q[idx_x[N][8]] = γ*(X[N][8] - ϵ_f)
    # rr = normalize(xf[1:3])
    # Qn = 1000*(I - rr*rr')
    # P[idx_x[N][1:3],idx_x[N][1:3]] = Qn'*Qn
    # q[idx_x[N][1:3]] = -(Qn'*Qn)'*(xf[1:3] - X[N][1:3])

    # solve using OSQP
    Q = copy(P)
    q = copy(q)
    A = copy(A_eq)
    b = copy(up_eq)
    G = [A_ineq;-A_ineq]
    h = [up_tr;-low_tr]

    z = quadprog(Q,q,A,b,G,h)
    δx = [z[idx_x[i]] for i = 1:(N)]
    δu = [z[idx_u[i]] for i = 1:(N-1)]

    # osqp = OSQP.Model()
    # OSQP.setup!(osqp; P=P, q=q, A=A, l=(Lo), u=(Up), eps_abs = 1e-8, eps_rel = 1e-8, max_iter = 10000,polish = 1)
    # results = OSQP.solve!(osqp)

    # # recover solution
    # δx = [results.x[idx_x[i]] for i = 1:(N)]
    # δu = [results.x[idx_u[i]] for i = 1:(N-1)]

    # apply the δ's
    cX = X + δx
    cU = U + δu

    return cU
end
