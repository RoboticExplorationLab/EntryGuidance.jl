using Convex, Mosek, MosekTools
using Attitude, LinearAlgebra, ForwardDiff, MATLAB
using StaticArrays
using SparseArrays
struct EntryVehicle{T} <: TO.AbstractModel
    evmodel::EG.CartesianModel{T}
end


function evdynamics(model::EntryVehicle, x, u)
    #unpack control
    α = deg2rad(18)
    #unpack state
    r = x[1:3]
    v = x[4:6]
    V = norm(v)

    #atmospheric density
    atmo_uncert = 1
    ρ = atmo_uncert*atmospheric_density([r[1]+model.evmodel.planet.R,0,0], model.evmodel)

    #Calculate drag acceleration
    Cd = drag_coefficient(α, model.evmodel)
    A = model.evmodel.vehicle.A
    m = model.evmodel.vehicle.m
    D = 0.5*Cd*ρ*A*V*V/m

    #Calculate lift acceleration
    Cl = lift_coefficient(α, model.evmodel)
    L = 0.5*Cl*ρ*A*V*V/m

    #get gravity
    g = gravitational_acceleration([r[1]+model.evmodel.planet.R,0,0], model.evmodel)
    # @show r
    # @show g
    @assert g[2] ≈ 0.0
    @assert g[3] ≈ 0.0
    @assert g[1] < 0.0


    #Aerodynamic acceleration
    e1 = cross([1;0;0],v)
    e1 = e1/norm(e1)
    e2 = cross(v,e1)
    e2 = e2/norm(e2)
    a = -(D/norm(v))*v + L*u[1]*e1 + L*u[2]*e2

    v̇ = a + g #- 2*Ω̂*v - Ω̂*Ω̂*r

    return [v; v̇]
end
function rk4(model,x_n,u,dt)
    dt = u[3]/3600
    k1 = dt*evdynamics(model,x_n,u)
    k2 = dt*evdynamics(model,x_n+k1/2,u)
    k3 = dt*evdynamics(model,x_n+k2/2,u)
    k4 = dt*evdynamics(model,x_n+k3,u)
    return (x_n + (1/6)*(k1 + 2*k2 + 2*k3 + k4))
end
struct EntryVehicle{T} <: TO.AbstractModel
    evmodel::EG.CartesianModel{T}
end
model = EntryVehicle(CartesianMSLModel())

#Initial conditions for MSL
Rm = model.evmodel.planet.R
r0 = [125.0, 0.0, 0.0] #Atmospheric interface at 125 km altitude
V0 = 5.845*3600 #Mars-relative velocity at interface of 5.845 km/sec
γ0 = -15.474*(pi/180.0) #Flight path angle at interface
v0 = V0*[sin(γ0), cos(γ0), 0.0]
x0 = [r0; v0]

xf = [10,450, 8.0,0,0,0]
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

function santize(X,U,A,B)
    @assert (length(A) == length(B))
    @assert (length(X) == length(U)+1)
    @assert (length(A) == length(U))
end
function mpc(model::EntryVehicle,X,U,xf,dt)

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

    # goal constraint
    # push!(cons, X[N][1:3] + δx[1:3,end] == xf[1:3])


    # dynamics constraints
    for i = 1:N-1
        dynamics_res = X[i+1] - rk4(model,X[i],U[i],dt)
        push!(cons, δx[:,i+1] == A[i]*δx[:,i] + B[i]*δu[:,i] - dynamics_res)
    end

    # control constraints
    for i = 1:N-1
        push!(cons, norm(δu[1:2,i],2) <= 1.0)
        push!(cons,δu[3,i]>= 0.5)
        push!(cons,δu[3,i]<= 3.0)
    end

    p = 0
    for i = 1:N-1
        # regularize U
        # p += sumsquares(U[i][1:2] + δu[1:2,i])
        # p += sumsquares( δu[1:2,i])

        # cost the dt to encourage it to be around 2 seconds
        # p += sumsquares( (U[i][3] + δu[3,i]) - 2  )
    end

    p += 100*sumsquares(X[N][2:3] + δx[2:3,end] - xf[2:3])
    p += 100000*sumsquares(vec(δu))

    push!(cons, X[N][1] + δx[1,end] == 10  )
    # trust region stuff
    # push!(cons, norm(vec(δx))<=5.0)

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


    N = 80
    dt = 69

    # initial guess
    N = 80
    X = [zeros(6) for i = 1:N]
    U = [[0;0.8;2] for i =1:N-1]
    X[1] = deepcopy(x0)

    for i = 1:(N-1)
        X[i+1] = rk4(model,X[i],U[i],dt)
    end

    cc = c_val(model,X,U,dt)

    TT = 25
    c_hist = zeros(TT)
    dx_hist = zeros(TT)
    @show norm(cc)
    # error()

    for i = 1:TT
        X2, U2 = mpc(model::EntryVehicle,X,U,xf,dt)

        @show norm(X - X2)
        dx_hist[i] = norm(X-X2)
        X = deepcopy(X2)
        U = deepcopy(U2)

        cc =  norm(c_val(model,X,U,dt))
        c_hist[i] = cc
        # dx_hist[i] = norm(X-X2)
        # @infiltrate
        # error()
        @show cc
        if cc < 1e-4
            @info "success"
            break
        end
    end

    # @infiltrate
    # error()
    mat"
    figure
    hold on
    plot($c_hist)
    hold off
    "
    # @infiltrate
    mat"
    figure
    hold on
    plot($dx_hist)
    hold off
    "
    # Xs = [zeros(6) for i = 1:N]
    # Xs[1] = x0
    # # rollout to see how solution looks
    # for i = 1:(N-1)
    #     Xs[i+1] = rk4(model,Xs[i],U[i],dt)
    # end
    #
    # xm = mat_from_vec(Xs)
    # um = mat_from_vec(U)
    # mat"
    # figure
    # hold on
    # plot($xm')
    # hold off
    # "
    # mat"
    # figure
    # hold on
    # stairs($um')
    # hold off
    # "



end

runsqp()
