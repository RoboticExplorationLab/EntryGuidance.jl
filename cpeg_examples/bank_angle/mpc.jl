using Convex, Mosek, MosekTools


function eg_mpc_quad(model::EntryVehicle,A,B,X,U,xf)

    # sizes for state and control
    nx = 7
    nu = 1
    N = length(X)

    # indicees for state and control
    idx_x = [(i-1)*(nx+nu) .+ (1:nx) for i = 1:length(X)]
    idx_u = [(i-1)*(nx+nu) .+ nx .+ (1:nu) for i = 1:(length(X)-1)]

    # constraint jacobian (in this form L≦Az≤U)
    nz = (N*nx) + (N-1)*nu
    nc = (N-1)*nx + nx + (N)*1
    A_eq = spzeros(nc,nz)

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
        A_eq[N*nx + i,idx_x[i][7]] = 1
    end

    # constraint bounds
    low_eq = zeros((N-1)*nx)
    up_eq = copy(low_eq)
    low_x0 = zeros(nx)
    up_x0 = zeros(nx)
    up_tr = deg2rad(20)*ones(N)
    low_tr = -deg2rad(20)*ones(N)

    # stack everything up
    A = A_eq
    Lo = [low_eq;low_x0;low_tr]
    Up = [up_eq;up_x0;up_tr]

    # cost function terms
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

    return cX, cU, norm(mat_from_vec(δu))
end

function eg_mpc_l1(model::EntryVehicle,A,B,X,U,xf)

    # dimension checking (debugging)
    @assert (length(A) == length(B))
    @assert (length(X) == length(U)+1)
    @assert (length(A) == length(U))
    N = length(X)
    x0 = copy(X[1])

    # initial condition
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

    # trust region on change to bank angle
    for i = 1:N
        push!(cons, δx[7,i] <=  deg2rad(20))
        push!(cons, δx[7,i] >= -deg2rad(20))
    end

    # cost function
    β = 1/length(U)
    problem = minimize( 1000*sumsquares( ( Qn*(X[N][1:3] + δx[1:3,N]) - xf[1:3])  )  + β*norm(vec(mat_from_vec(U)) + δu,1) , cons)

    # solve with Mosek
    Convex.solve!(problem, () -> Mosek.Optimizer())

    # recover solution
    cX = vec_from_mat(mat_from_vec(X) + evaluate(δx))
    cU = deepcopy(U)
    for i = 1:length(U)
        cU[i][1] += δu.value[i][1]
    end

    return cX, cU, norm(evaluate(δu))
end
