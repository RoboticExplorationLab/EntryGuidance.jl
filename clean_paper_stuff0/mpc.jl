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


    α = 1
    rr = normalize(xf[1:3])
    Qn = I - rr*rr'

    problem = minimize( norm( Qn*( (X[N][1:3] + δx[1:3,N]) - xf[1:3])  ) +  α*sumsquares(vec(δu)) , cons)

    # be gentle with the first 5 control inputs

    # problem = minimize( norm( Qn*( (X[N][1:3] + δx[1:3,N]) - xf[1:3])  ) , cons)
    Convex.solve!(problem, () -> Mosek.Optimizer())



    cX = vec_from_mat(mat_from_vec(X) + evaluate(δx))
    cU = vec_from_mat(mat_from_vec(U) + evaluate(δu))

    # Xm = mat_from_vec(X) + evaluate(δx)
    # Um = mat_from_vec(U) + evaluate(δu)
    #
    # dxm = evaluate(δx)
    # dum = evaluate(δu)
    #
    # mat"
    # figure
    # hold on
    # title('X')
    # plot($Xm')
    # hold off
    # "
    #
    # mat"
    # figure
    # hold on
    # title('U')
    # plot($Um')
    # hold off
    # "
    #
    # mat"
    # figure
    # hold on
    # title('dX')
    # plot($dxm')
    # hold off
    # "
    #
    # mat"
    # figure
    # hold on
    # title('dU')
    # plot($dum')
    # hold on
    # "
    #
    # @infiltrate
    # error()

    return cX, cU

end
