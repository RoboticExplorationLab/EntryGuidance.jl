using Convex, Mosek, MosekTools



function eg_mpc(A,B,X,U)

    @assert (length(A) == length(B))
    @assert (length(X) == length(U)+1)
    @assert (length(A) == length(U))
    N = length(X)
    x0 = copy(X[1])

    xm = mat_from_vec(X)
    um = mat_from_vec(U)

    mpc = (A = A, B = B, X = X, U = U)
    # @save "mpc.jld2" mpc
    mat"
    figure
    hold on
    title('Pre-Positions')
    plot($xm(1:3,:)')
    hold off
    "
    mat"
    figure
    hold on
    title('Pre-Velocities')
    plot($xm(1:3,:)')
    hold off
    "
    mat"
    figure
    hold on
    title('Pre-controls')
    plot($um')
    hold off
    "
    # max lift vector for each time
    L_max = zeros(N)
    for i = 1:N

        # model stuff
        Area = 1.0  # m²
        m = 1000 # kg
        ρ = 1.2*10 # kg/m³
        # Cl = 1.42*α (in radians)
        Cl = 1.42*deg2rad(30)

        # velocity
        V = norm(X[i][4:6])

        # this is maximum allowable lift
        L_max[i] = 0.5*Cl*ρ*Area*V*V/m
    end

    # starting δx
    δx0 = x0 - X[1]
    @assert norm(x0 - X[1]) ≈ 0
    # variables
    δx = Variable(6,N)
    δu = Variable(2,N-1)

    # dynamics constraints
    cons = Constraint[ δx[:,1]==δx0 ]

    for i = 1:N-1
        push!(cons, δx[:,i+1] == sparse(A[i])*δx[:,i] + B[i]*δu[:,i])
    end

    # lift vector constraints
    for i = 1:N-1
        push!(cons, norm(U[i] + δu[:,i]) <= L_max[i])
        # push!(cons, sumsquares(U[i] + δu[:,i]) <= L_max[i]^2)
    end

    # problem = minimize(norm(X[N][1:2] + δx[1:2,N]), cons)
    α = 1e1
    β = 1e2
    # problem = minimize(α*norm(X[N][1:2] + δx[1:2,N]) + norm(vec(mat_from_vec(U) + δu)) + norm(vec(δu)), cons)
    problem = minimize(α*norm(X[N][1:2] + δx[1:2,N]) + norm(vec(mat_from_vec(U) + δu)), cons)
    # problem = minimize(α*sumsquares(X[N][1:2] + δx[1:2,N]) + sumsquares(vec(mat_from_vec(U) + δu)), cons)
    # Convex.solve!(problem, () -> Mosek.Optimizer())
    # Convex.solve!(problem, () -> ECOS.Optimizer(abstol = 1e-9))
    # Convex.solve!(problem, () -> Gurobi.Optimizer())
    # Convex.solve!(problem, () -> COSMO.Optimizer(verbose = false))
    Convex.solve!(problem, () -> COSMO.Optimizer(eps_abs = 1e-12))

    cX = vec_from_mat(mat_from_vec(X) + evaluate(δx))
    cU = vec_from_mat(mat_from_vec(U) + evaluate(δu))

    dx = δx.value
    du = δu.value

    naiveJ = α*norm(X[N][1:2]) + norm(vec(mat_from_vec(U)))
    theirJ = α*norm(X[N][1:2] + dx[1:2,N]) + norm(vec(mat_from_vec(U) + du))
    @info "naive J is $naiveJ"
    @info "their J is $theirJ"
    mat"
    figure
    hold on
    title('Delta Position')
    plot($dx(1:3,:)')
    hold off
    "
    mat"
    figure
    hold on
    title('Delta Velocity')
    plot($dx(4:6,:)')
    hold off
    "
    mat"
    figure
    hold on
    title('Delta Control')
    plot($du')
    hold off
    "
    @infiltrate
    @info δx.value[1:3,end]
    # @infiltrate
    # error()
    return cX, cU
end
