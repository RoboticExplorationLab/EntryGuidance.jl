using Convex, Mosek, MosekTools



function eg_mpc(model::EntryVehicle,A,B,X,U,xf)

    @assert (length(A) == length(B))
    @assert (length(X) == length(U)+1)
    @assert (length(A) == length(U))
    N = length(X)
    x0 = copy(X[1])

    # get maximum lift
    L_max = zeros(N-1)
    for i = 1:N-1
        L_max[i] = getmaxL(model,X[i])
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


    α = 1e-2
    β = 0e1
    rr = normalize(xf[1:3])
    Qn = I - rr*rr'
    l1val = 0
    # for i = 1:(N-2)
    #     l1val += norm( ((U[i+1] + δu[:,i+1])/L_max[i+1]) - ((U[i] + δu[:,i])/L_max[i]) )
    # end
    γ = 1e1
    # problem = minimize( norm( Qn*( (X[N][1:3] + δx[1:3,N]) - xf[1:3])  ) + β*l1val +  α*sumsquares(vec(δu)) + γ*norm(X[N][4:6]), cons)

    # trust region stuff
    # push!(cons, norm(δx[:,N])<=200.0)

    # g loading constraint
    # for i = 1:N-1
    #     v1 = (X[i][4:6] + δx[4:6,i])*1000/3600
    #     v2 = (X[i+1][4:6] + δx[4:6,i+1])*1000/3600
    #     push!(cons,norm(v2 - v1)/(2*9.8) <= 13.0)
    # end # dt = 2.0 s

    # dynamic pressure constraint
    # for i = 1:N
    #     gk = pressure_constraint(model,X[i])
    #     ∇gk = pressure_gradient(model,X[i])
    #     push!(cons, dot(∇gk,δx[:,i]) <= -gk)
    # end
    problem = minimize( norm( Qn*( (X[N][1:3] + δx[1:3,N]) - xf[1:3])  ) +  α*sumsquares(vec(δu)) , cons)
    # problem = minimize( norm( ( (X[N][1:3] + δx[1:3,N]) - xf[1:3])  ) + γ*norm( vec(   mat_from_vec(U) + δu       )), cons)
    Convex.solve!(problem, () -> Mosek.Optimizer())

    cX = vec_from_mat(mat_from_vec(X) + evaluate(δx))
    cU = vec_from_mat(mat_from_vec(U) + evaluate(δu))
    return cX, cU


    # # JuMP version
    # @assert (length(A) == length(B))
    # @assert (length(X) == length(U)+1)
    # @assert (length(A) == length(U))
    # N = length(X)
    # x0 = copy(X[1])
    #
    # jmodel = Model(Mosek.Optimizer)
    # # set_optimizer_attribute(jmodel, "MSK_DPAR_INTPNT_CO_TOL_DFEAS",1e-15)
    # # set_optimizer_attribute(jmodel, "MSK_DPAR_INTPNT_CO_TOL_REL_GAP",1e-15)
    #
    # # jmodel = Model(COSMO.Optimizer)
    # # set_optimizer_attribute(jmodel, "check_termination",30)
    #
    #
    # @variable(jmodel, δx[1:6,1:N])
    # @variable(jmodel, δu[1:2,1:N-1])
    #
    # @constraint(jmodel,δx[:,1] .== zeros(6))
    # for i = 1:N-1
    #     @constraint(jmodel, δx[:,i+1] .== sparse(A[i])*δx[:,i] + sparse(B[i])*δu[:,i])
    # end
    #
    # # get maximum lift
    # maxL = zeros(N-1)
    # for i = 1:N-1
    #     maxL[i] = getmaxL(model,X[i])
    # end
    #
    # for i = 1:N-1
    #     @constraint(jmodel, [maxL[i],U[i][1] + δu[1,i],U[i][2] + δu[2,i]] in SecondOrderCone())
    # end
    #
    # # rr = normalize(xf[1:3])
    # # Qn = I - rr*rr'
    # # Q = Qn'*Qn
    # Q = I(3)
    # objective_exp = @expression(jmodel, (X[N][1:3] + δx[1:3,N] - xf[1:3] )'*(X[N][1:3] + δx[1:3,N] - xf[1:3] ))
    #
    # α = 1
    # for i = 1:N-1
    #     add_to_expression!(objective_exp, α*δu[:,i]'*δu[:,i])
    # end
    #
    # @objective(jmodel, Min, objective_exp)
    #
    # optimize!(jmodel)
    #
    # du = value.(δu)
    # dx = value.(δx)
    #
    # cX = vec_from_mat(mat_from_vec(X) + dx)
    # cU = vec_from_mat(mat_from_vec(U) + du)
    #
    # return cX, cU
end
