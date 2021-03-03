using Interpolations

struct EntryVehicle_fixed_time{T} <: TO.AbstractModel
    evmodel::EG.CartesianModel{T}
end
function RD.dynamics(model::EntryVehicle_fixed_time, x, u)
    α = u[2]
    σ̇ = u[1]
    return [EG.dynamics(x[1:6], EG.angles_input([α, x[7]],x[1:6],model.evmodel), model.evmodel); σ̇]
end

# Define a custom integration method
abstract type EntryVehicleRK_fixed_time <: RD.Explicit end

# Define the discrete dynamics function
function RD.discrete_dynamics(::Type{EntryVehicleRK_fixed_time}, model::EntryVehicle_fixed_time,
        x::StaticVector, u::StaticVector, t, dt)

    # h = u[3]/3600.0 #u is in seconds, dynamics are in hours
    h = dt/3600 # dt is seconts, h is in hours

    k1 = RD.dynamics(model, x,             u)*h;
    k2 = RD.dynamics(model, x + k1/2,      u)*h;
    k3 = RD.dynamics(model, x - k1 + 2*k2, u)*h;

    return x + (k1 + 4*k2 + k3)/6
end

Base.size(::EntryVehicle_fixed_time) = 7,2

function getK(X,U,t_traj,Q,R)

    x_int = LinearInterpolation(t_traj,X)
    u_int = LinearInterpolation(t_traj[1:end-1],U)

    dt = 1.0
    new_t = 0:dt:250

    N = length(new_t)

    new_x = x_int(new_t)
    new_u = u_int(new_t)



    model = EntryVehicle_fixed_time(CartesianMSLModel())
    n,m = size(model)

    A_vec = [@SMatrix zeros(n,n) for i = 1:N]
    B_vec = [@SMatrix zeros(n,m) for i = 1:N]

    AB = zeros(n,n+m)
    for i = 1:N
        x = new_x[i]
        u = new_u[i]
        z = KnotPoint(x,u[1:2],dt)
        discrete_jacobian!(EntryVehicleRK_fixed_time,AB,model,z)
        A_vec[i] = SMatrix{n,n}(AB[:,1:n])
        B_vec[i] = SMatrix{n,m}(AB[:,(n+1):end])
    end

    K = [@SMatrix zeros(m,n) for i = 1:N]
    S = [@SMatrix zeros(n,n) for i = 1:N]

    S[N] = Q
    for k = (N-1):(-1):1
        Ak = A_vec[k]
        Bk = B_vec[k]
        K[k] = (R+Bk'*S[k+1]*Bk)\(Bk'*S[k+1]*Ak)
        S[k] = Q + K[k]'*R*K[k] + (Ak - Bk*K[k])'*S[k+1]*(Ak - Bk*K[k])
    end

    # get a K object
    K_itp = interpolate(K, BSpline(Linear())) #works
    K_int = scale(K_itp, new_t)

    return new_t, x_int, u_int, K_itp, K
end
