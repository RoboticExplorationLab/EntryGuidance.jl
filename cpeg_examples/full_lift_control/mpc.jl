using Convex, Mosek, MosekTools


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

    γ = 10000         # miss distance penalty
    α = 1e-6/length(U) # regularizer
    β = 1e-6/length(U) # control penalty


    problem = minimize(γ*sumsquares( Qn*( (X[N][1:3] + δx[1:3,N]) - xf[1:3])  ) +  α*sumsquares(vec(δu)) + β*sumsquares(vec(mat_from_vec(U)) + vec(δu)), cons)

    # solve with Mosek
    Convex.solve!(problem, () -> Mosek.Optimizer())

    # recover solutions
    cX = vec_from_mat(mat_from_vec(X) + evaluate(δx))
    cU = vec_from_mat(mat_from_vec(U) + evaluate(δu))

    return cX, cU, norm(evaluate(δu))

end
