using Convex, Mosek, MosekTools
using Attitude, LinearAlgebra, ForwardDiff, MATLAB
using StaticArrays

function rk4(model,x_n,u,dt)
    k1 = dt*dynamics(model,x_n,u)
    k2 = dt*dynamics(model,x_n+k1/2,u)
    k3 = dt*dynamics(model,x_n+k2/2,u)
    k4 = dt*dynamics(model,x_n+k3,u)
    return (x_n + (1/6)*(k1 + 2*k2 + 2*k3 + k4))
end
function getAB(model,X,U,dt)
    N = length(X)
    n = length(X[1])
    m = length(U[1])
    A = [spzeros(n,n) for i = 1:(N-1)]
    B = [spzeros(n,m) for i = 1:(N-1)]
    for k = 1:(N-1)
        A[k] = sparse(ForwardDiff.jacobian(_x -> rk4(model,_x,U[k],dt),X[k]))
        B[k] = sparse(ForwardDiff.jacobian(_u -> rk4(model,X[k],_u,dt),U[k]))
    end
    santize(X,U,A,B)
    return A,B
end

struct SAT
    J::Diagonal{Float64,SArray{Tuple{3},Float64,1,3}}
end

function dynamics(model::SAT,x,u)
    p = x[SA[1,2,3]]
    ω = x[SA[4,5,6]]
    α = model.J\(u - cross(ω,model.J*ω))
    pdot = pdot_from_w(p,ω)
    return [pdot[1],pdot[2],pdot[3],α[1],α[2],α[3]]
end
function santize(X,U,A,B)
    @assert (length(A) == length(B))
    @assert (length(X) == length(U)+1)
    @assert (length(A) == length(U))
end
function mpc(model::SAT,X,U,xf,dt)

    # linearize the system
    A,B = getAB(model,X,U,dt)

    # problem length and IC
    N = length(X)
    x0 = copy(X[1])

    # variables
    δx = Variable(6,N)
    δu = Variable(3,N-1)

    # initial condition of all zeros for δx
    cons = Constraint[ δx[:,1]==zeros(length(X[1])) ]

    # dynamics constraints
    for i = 1:N-1
        dynamics_res = X[i+1] - rk4(model,X[i],U[i],dt)
        push!(cons, δx[:,i+1] == A[i]*δx[:,i] + B[i]*δu[:,i] - dynamics_res)

        push!(cons, U[i] + δu[:,i] <=  ones(3))
        push!(cons, U[i] + δu[:,i] >= -ones(3))
    end

    # trust region stuff
    push!(cons, norm(vec(δx))<=2.0)

    # LQR cost
    Q = Diagonal(ones(6))
    R = Diagonal(ones(3))
    p = 0
    # for i = 1:N-1
    #     p += quadform((X[i] + δx[:,i] - xf),Q)
    #     p += quadform(U[i] + δu[:,i],R)
    # end
    # p += quadform((X[N] + δx[:,N] - xf),Q)
    for i = 1:N-1
        p += norm((X[i] + δx[:,i] - xf),2)
        p += norm(U[i] + δu[:,i],1)
    end
    p += norm((X[N] + δx[:,N] - xf))

    # solve problem
    problem = minimize(p, cons)
    Convex.solve!(problem, () -> Mosek.Optimizer())
    cX = vec_from_mat(mat_from_vec(X) + evaluate(δx))
    cU = vec_from_mat(mat_from_vec(U) + evaluate(δu))
    return cX, cU
end

function c_val(model,X,U,dt)
    c = zeros((length(X)-1)*length(X[1]))
    nx = length(X[1])
    N = length(X)
    for i = 1:(length(X)-1)
        c[ nx*(i-1) .+ (1:nx)  ] = X[i+1] - rk4(model,X[i],U[i],dt)
    end
    return c
end
function runsqp()


    N = 50
    dt = .5

    x0 = [p_from_phi(deg2rad(170)*normalize([1;2;3])); zeros(3)]
    xf = zeros(6)
    X = [copy(x0) for i = 1:N]
    U = [zeros(3) for i = 1:(N-1)]

    model = SAT(Diagonal(SA[1,2,3.0]))

    for i = 1:10
        X2, U2 = mpc(model::SAT,X,U,xf,dt)

        @show norm(X - X2)
        X = deepcopy(X2)
        U = deepcopy(U2)

        @info norm(c_val(model,X,U,dt))
    end

    Xs = [zeros(6) for i = 1:N]
    Xs[1] = x0
    # rollout to see how solution looks
    for i = 1:(N-1)
        Xs[i+1] = rk4(model,Xs[i],U[i],dt)
    end

    xm = mat_from_vec(Xs)
    um = mat_from_vec(U)
    mat"
    figure
    hold on
    plot($xm')
    hold off
    "
    mat"
    figure
    hold on
    stairs($um')
    hold off
    "



end

runsqp()
