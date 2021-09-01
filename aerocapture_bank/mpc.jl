using Convex, Mosek, MosekTools
using JuMP


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
    # push!(cons, norm(δx[7,:],Inf) <= deg2rad(10))
    for i = 1:N
        push!(cons, abs(δx[7,:]) <= deg2rad(10))
    end

    problem = minimize(γ*p +  α*sumsquares(vec(δu)), cons)

    Convex.solve!(problem, () -> Mosek.Optimizer())

    cX = vec_from_mat(mat_from_vec(X) + evaluate(δx))
    cU = vec_from_mat(mat_from_vec(U) + evaluate(δu))

    # @infiltrate
    return cU

end
