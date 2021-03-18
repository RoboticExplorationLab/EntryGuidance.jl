using Convex, Mosek, MosekTools



function eg_mpc(Af,Bf,Xf,Uf,x0,idx,N,model)

    # trim to be the correct indices
    A = Af[idx:(idx+N)]
    B = Bf[idx:(idx+N)]
    X = Xf[idx:(idx+N)]
    U = Uf[idx:(idx+N)]

    # vehicle stuff
    A = model.evmodel.vehicle.A
    m = model.evmodel.vehicle.m

    # max lift vector for each time
    L = [zeros(2) for i = 1:N]
    L_max = zeros(N)
    for i = 1:N
        Cl_max = .619 # 25 degrees
        ρ = amtospheric_density(X[i],model.evmodel)
        V = norm(X[i][4:6])
        L_max = 0.5*Cl_max*ρ*A*V*V/m

        # actual control inputs
        α,σ = U[i]
        Cl = lift_coefficient(α, model.evmodel)
        L_mag = 0.5*Cl*ρ*A*V*V/m
        L[i] = L_mag*[sin(σ); cos(σ)]
    end

    # starting δx
    δx0 = x0 - X[1]

    # variables
    δx = Variable(6,N)
    δu = Variable(2,N)

    # dynamics constraints
    cons = Constraint[ δx[:,1]==δx0 ]
    for i = 1:N-1
        push!(cons, δx[i+1] == A[i]*δx[i] + B[i]*δu[i])
    end

    # lift vector constraints
    for i = 1:N
        push!(norm(U[i] + δu[i]) <= L_max[i])
    end

    problem = minimize(sumsquares(δx), constraints)
