using LinearAlgebra
using StaticArrays
using EntryGuidance
const EG = EntryGuidance
using TrajectoryOptimization
const TO = TrajectoryOptimization
using RobotDynamics
const RD = RobotDynamics
using Altro
using Test
using Attitude

struct Sat2 <: TO.AbstractModel
    J::Diagonal{Int64,SArray{Tuple{3},Int64,1,3}}
    invJ::Diagonal{Float64,SArray{Tuple{3},Float64,1,3}}
end
function RD.dynamics(model::Sat2, x, u)
    p = SA[x[1],x[2],x[3]]
    ω = SA[x[4],x[5],x[6]]
    pdot = pdot_from_w(p,ω)
    α = model.invJ*(u - cross(ω,model.J*ω))
    return SA[pdot[1],pdot[2],pdot[3],α[1],α[2],α[3]]
end


Base.size(::Sat2) = 6,3

J = Diagonal(@SVector [1,2,3])::Diagonal{Int64,SArray{Tuple{3},Int64,1,3}}
invJ = inv(J)::Diagonal{Float64,SArray{Tuple{3},Float64,1,3}}
Sat2() = Sat2(J,invJ)
model = Sat2()
n,m = size(model)
tf = 30.0 #time is in hours
N = 121
dt0 = (tf)/(N-1)

#Initial conditions for MSL
x0 = SVector{n}([p_from_phi(deg2rad(110)*normalize(randn(3))); zeros(3)])
xf = @SVector zeros(n)

#Cost function
Q = Diagonal(SVector{n}(ones(n)))
R = Diagonal(SVector{m}(ones(m)))

stage_cost = LQRCost(Q,R,xf,terminal=false)

Qn = 10*Q
terminal_cost = LQRCost(Qn,R,xf,terminal=true)

obj = Objective(stage_cost,terminal_cost,N)

#Constraints
cons = TO.ConstraintList(n,m,N)
# ∞ = Inf64
# add_constraint!(cons, BoundConstraint(n,m,x_min=SA[-∞,-∞,-∞,-∞,-∞,-∞,-pi],x_max=SA[∞,∞,∞,∞,∞,∞,pi]),1:N)
add_constraint!(cons, BoundConstraint(n,m,u_min=-.05*SA[1,1,1],u_max=.05*SA[1,1,1]),1:(N-1))
# add_constraint!(cons, GoalConstraint(xf, [1,2,3]), N:N)

#Initial Controls
u_traj = ones(m,N-1)
u_traj[2,:] .= 15.0*pi/180
u_traj[3,:] .= dt0*ones(N-1)

prob = TO.Problem(model, obj, xf, tf, x0=x0, constraints=cons, integration=RK4)

solver = ALTROSolver(prob)
# solver.opts.max_cost_value = 1e15
# solver.opts.bp_reg_initial = 1e-6
# solver.opts.bp_reg_min = 1e-6
# solver.opts.constraint_tolerance = 1e-2
# solver.opts.cost_tolerance_intermediate = 1e-4
solver.opts.projected_newton = false
solver.opts.verbose = 2
solve!(solver)
#
X = states(solver)
U = controls(solver)
#
using MATLAB

Xm = mat_from_vec(X)
Um = mat_from_vec(U)
mat"
figure
hold on
plot($Xm')
hold off
"

mat"
figure
hold on
plot($Um')
hold off
"

function TVLQR(model,X,U,dt,integration,Q,Qf,R)

    N = length(X)
    n,m = size(model)
    A = [zeros(n,n) for i = 1:N-1]
    B = [zeros(n,m) for i = 1:N-1]

    # get jacobians along the trajectory
    ∇f = zeros(n,n+m)
    for i = 1:N-1
        z = KnotPoint(X[i],U[i],dt)
        discrete_jacobian!(integration,∇f,model,z)
        A[i] = ∇f[1:n,1:n]
        B[i] = ∇f[1:n,n+1:end]
    end

    S = [zeros(n,n) for i = 1:N]
    K = [zeros(m,n) for i = 1:N-1]
    S[N] = Qf
    for k = (N-1):(-1):1
        # get views
        Bk = B[k]
        Ak = A[k]
        Skp1 = S[k+1]

        # feedback solve
        K[k] = (R + Bk'*Skp1*Bk)\Bk'*Skp1*Ak

        # cost to go
        Kk = K[k]
        S[k] = Q + Kk'*R*Kk + (Ak - Bk*Kk)'*Skp1*(Ak - Bk*Kk)
    end
    return K
end


Kt = TVLQR(model,X,U,dt0,RK4,Q,Qn,R/100)


Xs = [@SVector zeros(n) for i = 1:N]
Us = [@SVector zeros(m) for i = 1:N]
# tscale = 3600
dt = dt0
Xs[1] = x0

for i = 1:(N-1)
    dx = Xs[i]  - X[i]
    # Us[i] = U[i] - Kt[i]*dx
    Us[i] = U[i] + [.001;.002;.003]- Kt[i]*dx
    z = KnotPoint(Xs[i],Us[i],dt)
    Xs[i+1] = discrete_dynamics(RK4,model,z) #+ .005*randn(n)
end

Xsm = mat_from_vec(Xs)

mat"
figure
title('MRP')
hold on
plot($Xm(1:3,:)')
plot($Xsm(1:3,:)','o')
hold off
"
mat"
figure
title('Angular Velocity')
hold on
plot($Xm(4:6,:)')
plot($Xsm(4:6,:)','o')
hold off
"
# mat"
# figure
# hold on
# plot($Um')
# hold off
# "
