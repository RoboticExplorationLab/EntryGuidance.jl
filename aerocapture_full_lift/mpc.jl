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
    δx = Variable(7,N)
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

    γ = 1e5         # miss distance penalty
    α = 1 # regularizer
    β = 1 # control penalty

    # p = 0.0
    # for i = 1:N
    #     ϵ = X[i][7] + δx[7,i]
    #     p+= square(ϵ - ϵ_f)
    # end
    # Xm = mat_from_vec(X)
    p = square(X[N][7] + δx[7,N] - ϵ_f)


    # @infiltrate
    # error()
    # p = sumsquares( vec(Xm[7,:]) + vec(δx[7,:]) - ϵ_f*ones(length(X)) )
    # p = 1
    problem = minimize(γ*p +  α*sumsquares(vec(δu)) + β*sumsquares(vec(mat_from_vec(U)) + vec(δu)), cons)

    Convex.solve!(problem, () -> Mosek.Optimizer())

    cX = vec_from_mat(mat_from_vec(X) + evaluate(δx))
    cU = vec_from_mat(mat_from_vec(U) + evaluate(δu))

    # @infiltrate
    return cU

end
