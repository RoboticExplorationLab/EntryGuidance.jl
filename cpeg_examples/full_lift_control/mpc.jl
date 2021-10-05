using Convex, Mosek, MosekTools
# using JuMP


function eg_mpc(model::EntryVehicle,A,B,X,U,xf, mpc_iteration)

    # quick size cehcks
    @assert (length(A) == length(B))
    @assert (length(X) == length(U)+1)
    @assert (length(A) == length(U))
    N = length(X)
    x0 = copy(X[1])

    # starting δx
    δx0 = x0 - X[1]
    @assert norm(x0 - X[1]) ≈ 0
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
    # @variable(jmodel, γ)
    #
    # @constraint(jmodel,δx[:,1] .== δx0)
    # for i = 1:N-1
    #     @constraint(jmodel, δx[:,i+1] .== sparse(A[i])*δx[:,i] + sparse(B[i])*δu[:,i])
    # end
    #
    #
    # for i = 1:N-1
    #     @constraint(jmodel, [1.0, U[i][1] + δu[1,i],U[i][2] + δu[2,i]] in SecondOrderCone())
    # end
    #
    # # add slack variable such that
    # # norm(Qn*(X[N][1:3] + δx[1:3,N] - xf[1:3]) ≦ γ
    # # then we cost gamma
    # rr = normalize(xf[1:3])
    # Qn = I - rr*rr'
    # @constraint(jmodel, [γ; Qn*(X[N][1:3] + δx[1:3,N] - xf[1:3] )] in SecondOrderCone())
    #
    # # trust region
    # @constraint(jmodel, [5.0; Qn*δx[1:3,N]] in SecondOrderCone())
    #
    # α = 1
    # β = 1e-3
    # P = α*I(length(vec(δu)))
    #
    # Uvec = vec(mat_from_vec(U))
    # # @objective(jmodel, Min, γ + α*vec(δu)'*vec(δu) + β*(Uvec + vec(δu))'*(Uvec + vec(δu)))
    # @objective(jmodel, Min, γ + α*vec(δu)'*vec(δu))
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
        push!(cons, norm(U[i] + δu[:,i]) <= 1.0)
    end



    rr = normalize(xf[1:3])
    Qn = I - rr*rr'

    # trust region
    # push!(cons, norm( Qn*δx[1:3,N]) <= 5)
    # push!(cons,norm( Qn*( (X[N][1:3] + δx[1:3,N]) - xf[1:3])) <= norm( Qn*( (X[N][1:3]) - xf[1:3])  ) )

    γ = 100         # miss distance penalty
    α = 1/length(U) # regularizer
    β = 1/length(U) # control penalty


    problem = minimize(γ*norm( Qn*( (X[N][1:3] + δx[1:3,N]) - xf[1:3])  ) +  α*norm(vec(δu)) + β*sumsquares(vec(mat_from_vec(U)) + vec(δu)), cons)


    # problem = minimize( norm( Qn*( (X[N][1:3] + δx[1:3,N]) - xf[1:3])  ) , cons)
    Convex.solve!(problem, () -> Mosek.Optimizer())



    cX = vec_from_mat(mat_from_vec(X) + evaluate(δx))
    cU = vec_from_mat(mat_from_vec(U) + evaluate(δu))

    return cX, cU

end
