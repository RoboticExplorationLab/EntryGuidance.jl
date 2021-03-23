using Convex, Mosek, MosekTools



function eg_mpc(A,B,X,U,x0)

    @assert (length(A) == length(B))
    @assert (length(X) == length(U)+1)
    @assert (length(A) == length(U))
    N = length(X)

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

    # variables
    δx = Variable(6,N)
    δu = Variable(2,N-1)

    # dynamics constraints
    cons = Constraint[ δx[:,1]==δx0 ]
    # @infiltrate
    # error()
    for i = 1:N-1
        # @show i
        # @infiltrate
        # error()
        push!(cons, δx[:,i+1] == A[i]*δx[:,i] + B[i]*δu[:,i])
    end

    # lift vector constraints
    for i = 1:N-1
        push!(cons, norm(U[i] + δu[:,i]) <= L_max[i])
    end

    problem = minimize(norm(X[N][1:2] + δx[1:2,N]), cons)
    # @time Convex.solve!(problem, () -> Mosek.Optimizer())
    @time Convex.solve!(problem, () -> COSMO.Optimizer())

    cX = vec_from_mat(mat_from_vec(X) + evaluate(δx))
    cU = vec_from_mat(mat_from_vec(U) + evaluate(δu))
    # @infiltrate
    # error()
    return cX, cU
end
