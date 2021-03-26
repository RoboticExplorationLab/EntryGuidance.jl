using JuMP, Mosek, MosekTools


function jump_mpc(model::EntryVehicle,X,U,A,B,xf)

    @assert (length(A) == length(B))
    @assert (length(X) == length(U)+1)
    @assert (length(A) == length(U))
    N = length(X)
    x0 = copy(X[1])

    jmodel = Model(Mosek.Optimizer)
    set_optimizer_attribute(jmodel, "MSK_DPAR_INTPNT_CO_TOL_DFEAS",1e-15)
    set_optimizer_attribute(jmodel, "MSK_DPAR_INTPNT_CO_TOL_REL_GAP",1e-15)


    # get maximum lift
    maxL = zeros(N-1)
    for i = 1:N-1
        maxL[i] = getmaxL(model,X[i])
    end


    @variable(jmodel, δx[1:6,1:N])
    @variable(jmodel, δu[1:2,1:N-1])

    @constraint(jmodel,δx[:,1] .= zeros(6))
    for i = 1:N-1
        @constraint(jmodel, δx[:,i+1] .== A[i]*δx[:,i] + B[i]*δu[:,i])
    end
