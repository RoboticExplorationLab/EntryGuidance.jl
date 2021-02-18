using LinearAlgebra, RobotDynamics, Convex, Mosek, MosekTools
using StaticArrays
using Attitude
using Infiltrator

struct Sat <: AbstractModel
    J   ::Diagonal{Float64,SArray{Tuple{3},Float64,1,3}}
    invJ::Diagonal{Float64,SArray{Tuple{3},Float64,1,3}}
end

Sat() = Sat(Diagonal(SA[1,2,3.0]), inv(Diagonal(SA[1,2,3.0])))

function RobotDynamics.dynamics(model::Sat,x,u)
    J = model.J
    invJ = model.invJ

    p = x[SA[1,2,3]]
    ω = x[SA[4,5,6]]

    pdot = pdot_from_w(p,ω)
    ω̇ = invJ*(u - cross(ω,J*ω))

    return SA[pdot[1],pdot[2],pdot[3],ω̇[1],ω̇[2],ω̇[3]]
end

# Specify the state and control dimensions
RobotDynamics.state_dim(::Sat) = 6
RobotDynamics.control_dim(::Sat) = 3

# Create the model
model = Sat()
n,m = size(model)

# Generate random state and control vector
x,u = rand(model)
dt = 0.1  # time step (s)
z = KnotPoint(x,u,dt)

# # Evaluate the continuous dynamics and Jacobian
# ẋ = dynamics(model, x, u)
∇f = RobotDynamics.DynamicsJacobian(model)
# jacobian!(∇f, model, z)

# Evaluate the discrete dynamics and Jacobian
x′ = discrete_dynamics(RK3, model, z)
discrete_jacobian!(RK3, ∇f, model, z)

function jacobian(x,u,dt)
    z = KnotPoint(x,u,dt)
    n,m = size(model)
    H = zeros(n,n+m)
    discrete_jacobian!(RK3, H, model, z)
    return H[:,1:n], H[:,n+1:end]
end

# AAd, BBd = jacobian(x,u,dt)

N = 20
xtraj = [@SVector randn(n) for i = 1:N ]
utraj = [@SVector randn(m) for i = 1:N-1 ]

idx_c = [(i-1)*n .+ (1:n) for i = 1:(N-1)]
function dynacon(model::Sat, z_traj,dt,idx_x,idx_u,idx_c,N)
    c = zeros(eltype(z_traj),(N-1)*length(xtraj[1]))
    for i = 1:N-1

        # current step
        xk = z_traj[idx_x[i]]
        uk = z_traj[idx_u[i]]
        zk = KnotPoint(xk,uk,dt)

        # next state value
        xkp1 = z_traj[idx_x[i+1]]

        # x2 - f(x1,u1)
        c[idx_c[i]] = xkp1 - discrete_dynamics(RK3, model, zk)
    end
    return c
end

ztraj2 = [vec(mat_from_vec(xtraj));vec(mat_from_vec(utraj))]
c1 = dynacon(model::Sat, ztraj2,dt,idx_x,idx_u,idx_c,N)

idx_x = [(i-1)*n .+ (1:n) for i = 1:N]
idx_u = [( N*n + (i-1)*m) .+ (1:m) for i = 1:N-1]

using FiniteDiff
const FD2 = FiniteDiff
using ForwardDiff
const FD = ForwardDiff
_c(_z) = dynacon(model::Sat, _z,dt,idx_x,idx_u,idx_c,N)

# J = FD2.finite_difference_jacobian(_c,ztraj2)
J = FD.jacobian(_c,ztraj2)

function myJ(model::Sat, z_traj,dt,idx_x,idx_u,idx_c,N)
    J = zeros((N-1)*n,(N*n + (N-1)*m))

    for i = 1:(N-1)
        xk = z_traj[idx_x[i]]
        uk = z_traj[idx_u[i]]
        A,B = jacobian(xk,uk,dt)
        J[idx_c[i],idx_x[i]]   = -A
        J[idx_c[i],idx_u[i]]   = -B
        J[idx_c[i],idx_x[i+1]] = I(length(xk))
    end
    return J
end

J2 =  myJ(model::Sat, ztraj2,dt,idx_x,idx_u,idx_c,N)