using Convex, Mosek, MosekTools



function eg_mpc(model::EntryVehicle,A,B,X,U,xf)
    #
    # @assert (length(A) == length(B))
    # @assert (length(X) == length(U)+1)
    # @assert (length(A) == length(U))
    # N = length(X)
    # x0 = copy(X[1])
    #
    # # get maximum lift
    # L_max = zeros(N-1)
    # for i = 1:N-1
    #     L_max[i] = getmaxL(model,X[i])
    # end
    #
    # # starting δx
    # δx0 = x0 - X[1]
    # @assert norm(x0 - X[1]) ≈ 0
    # # variables
    # δx = Variable(6,N)
    # δu = Variable(2,N-1)
    #
    # # dynamics constraints
    # cons = Constraint[ δx[:,1]==δx0 ]
    #
    # for i = 1:N-1
    #     push!(cons, δx[:,i+1] == sparse(A[i])*δx[:,i] + B[i]*δu[:,i])
    # end
    #
    # # lift vector constraints
    # for i = 1:N-1
    #     push!(cons, norm(U[i] + δu[:,i]) <= L_max[i])
    #     # push!(cons, sumsquares(U[i] + δu[:,i]) <= L_max[i]^2)
    # end
    #
    #
    # α = 1e-1
    # β = 1e2
    # problem = minimize(α*norm(  (X[N][1:3] + δx[1:3,N]) - xf[1:3]  ) + norm(vec(δu)), cons)
    # # problem = minimize(α*norm(X[N][1:2] + δx[1:2,N]) + norm(vec(mat_from_vec(U) + δu)) + norm(vec(δu)), cons)
    # # problem = minimize(α*norm(X[N][1:2] + δx[1:2,N]) + norm(vec(mat_from_vec(U) + δu)), cons)
    # # problem = minimize(α*sumsquares(X[N][1:2] + δx[1:2,N]) + sumsquares(vec(mat_from_vec(U) + δu)), cons)
    # Convex.solve!(problem, () -> Mosek.Optimizer())
    # # Convex.solve!(problem, () -> ECOS.Optimizer(abstol = 1e-9))
    # # Convex.solve!(problem, () -> Gurobi.Optimizer())
    # # Convex.solve!(problem, () -> COSMO.Optimizer(verbose = false))
    # # Convex.solve!(problem, () -> COSMO.Optimizer(eps_abs = 1e-12))
    #
    # cX = vec_from_mat(mat_from_vec(X) + evaluate(δx))
    # cU = vec_from_mat(mat_from_vec(U) + evaluate(δu))
    #
    # dx = δx.value
    # du = δu.value
    # return cX, cU
    @assert (length(A) == length(B))
    @assert (length(X) == length(U)+1)
    @assert (length(A) == length(U))
    N = length(X)
    x0 = copy(X[1])

    jmodel = Model(Mosek.Optimizer)
    # set_optimizer_attribute(jmodel, "MSK_DPAR_INTPNT_CO_TOL_DFEAS",1e-15)
    # set_optimizer_attribute(jmodel, "MSK_DPAR_INTPNT_CO_TOL_REL_GAP",1e-15)

    # jmodel = Model(COSMO.Optimizer)
    # set_optimizer_attribute(jmodel, "check_termination",30)


    @variable(jmodel, δx[1:6,1:N])
    @variable(jmodel, δu[1:2,1:N-1])

    @constraint(jmodel,δx[:,1] .== zeros(6))
    for i = 1:N-1
        @constraint(jmodel, δx[:,i+1] .== A[i]*δx[:,i] + B[i]*δu[:,i])
    end

    # get maximum lift
    maxL = zeros(N-1)
    for i = 1:N-1
        maxL[i] = getmaxL(model,X[i])
    end

    for i = 1:N-1
        @constraint(jmodel, [maxL[i],δu[1,i],δu[2,i]] in SecondOrderCone())
    end

    # rr = normalize(xf[1:3])
    # Qn = I - rr*rr'
    # Q = Qn'*Qn
    Q = I(3)
    objective_exp = @expression(jmodel, (X[N][1:3] + δx[1:3,N] - xf[1:3] )' * Q * (X[N][1:3] + δx[1:3,N] - xf[1:3] ))

    α = 1e-2
    for i = 1:N-1
        add_to_expression!(objective_exp, α*δu[:,i]'*δu[:,i])
    end

    @objective(jmodel, Min, objective_exp)

    optimize!(jmodel)

    du = value.(δu)
    dx = value.(δx)

    cX = vec_from_mat(mat_from_vec(X) + dx)
    cU = vec_from_mat(mat_from_vec(U) + du)
    # @infiltrate
    # error()
    # return (X + vec_from_mat(dx)), (mat_from_vec(U) + du)
    return cX, cU
end
