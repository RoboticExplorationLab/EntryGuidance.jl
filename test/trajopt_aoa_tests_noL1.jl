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

struct EntryVehicle{T} <: TO.AbstractModel
    evmodel::EG.CartesianModel{T}
end

function RD.dynamics(model::EntryVehicle, x, u)
    α = u[2]
    σ̇ = u[1]
    return [EG.dynamics(x[1:6], EG.angles_input([α, x[7]],x[1:6],model.evmodel), model.evmodel); σ̇]
end

# Define a custom integration method
abstract type EntryVehicleRK <: RD.Explicit end

# Define the discrete dynamics function
function RD.discrete_dynamics(::Type{EntryVehicleRK}, model::EntryVehicle,
        x::StaticVector, u::StaticVector, t, dt)

    h = u[3]/3600.0 #u is in seconds, dynamics are in hours

    k1 = RD.dynamics(model, x,             u)*h;
    k2 = RD.dynamics(model, x + k1/2,      u)*h;
    k3 = RD.dynamics(model, x - k1 + 2*k2, u)*h;

    return x + (k1 + 4*k2 + k3)/6
end

Base.size(::EntryVehicle) = 7,3

model = EntryVehicle(CartesianMSLModel())
n,m = size(model)
tf = 5.5/60.0 #time is in hours
N = 121
dt0 = (tf*3600)/(N-1)

#Initial conditions for MSL
Rm = model.evmodel.planet.R
r0 = [Rm+125.0, 0.0, 0.0] #Atmospheric interface at 125 km altitude
V0 = 5.845*3600 #Mars-relative velocity at interface of 5.845 km/sec
γ0 = -15.474*(pi/180.0) #Flight path angle at interface
v0 = V0*[sin(γ0), cos(γ0), 0.0]
α0 = 15.0*pi/180 #derived from MSL initial L/D = 0.24
σ = 0.0*pi/180 #this is totally made up
x0 = SVector{n}([r0; v0; σ])

#Final conditions for MSL 631.979 km down range and 7.869 km cross-range from entry point
Rf = Rm+10.0 #Parachute opens at 10 km altitude
rf = Rf*[cos(7.869/Rf)*cos(631.979/Rf); cos(7.869/Rf)*sin(631.979/Rf); sin(7.869/Rf)]
vf = zeros(3)
xf = SVector{n}([rf; vf; σ])

#Cost function
Q = Diagonal(SVector{n}(zeros(n)))
R = Diagonal(SVector{m}([0.00001, 1, 0.1]))

stage_cost = LQRCost(Q,R,xf,terminal=false)

Qn = Diagonal(@SVector [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
terminal_cost = LQRCost(Qn,R,xf,terminal=true)

obj = Objective(stage_cost,terminal_cost,N)

#Constraints
cons = TO.ConstraintList(n,m,N)
∞ = Inf64
add_constraint!(cons, BoundConstraint(n,m,x_min=SA[-∞,-∞,-∞,-∞,-∞,-∞,-pi],x_max=SA[∞,∞,∞,∞,∞,∞,pi]),1:N)
add_constraint!(cons, BoundConstraint(n,m,u_min=[-∞,10.0*pi/180,2.0],u_max=[∞,20.0*pi/180,3.0]),1:(N-1))
add_constraint!(cons, GoalConstraint(xf, [1,2,3]), N:N)

#Initial Controls
u_traj = ones(m,N-1)
u_traj[2,:] .= 15.0*pi/180
u_traj[3,:] .= dt0*ones(N-1)

prob = TO.Problem(model, obj, xf, tf, x0=x0, U0=u_traj, constraints=cons, integration=EntryVehicleRK)

solver = ALTROSolver(prob)
solver.opts.max_cost_value = 1e15
solver.opts.bp_reg_initial = 1e-6
solver.opts.bp_reg_min = 1e-6
solver.opts.constraint_tolerance = 1e-2
solver.opts.cost_tolerance_intermediate = 1e-4
solver.opts.projected_newton = false
solver.opts.verbose = 2
solve!(solver)

X = states(solver)
U = controls(solver)

@test max_violation(solver) < 1e-2

# #Test rollout
x_traj = zeros(n,N)
# # x_traj[:,1] .= x0
# # for k = 1:(N-1)
# #     x_traj[:,k+1] .= discrete_dynamics(EntryVehicleRK,model,SVector{n}(x_traj[:,k]),SVector{m}(U0[:,k]),0.0,2.0)
# # end
alt = zeros(N)
bank = zeros(N)
#
for k = 1:N
    x_traj[:,k] .= X[k]
    alt[k] = norm(x_traj[1:3,k])-Rm
    bank[k] = x_traj[7,k]
end

u_traj = zeros(m,N-1)
AoA = zeros(N-1)
# # u_traj = U0
σ̇ = zeros(N-1)
dt = zeros(N-1)
t_traj = zeros(N)
down_range = zeros(N)
cross_range = zeros(N)
for k = 1:(N-1)
    u_traj[:,k] .= U[k]
    σ̇[k] = u_traj[1,k]
    AoA[k] = u_traj[2,k]
    dt[k] = u_traj[3,k]
    t_traj[k+1] = t_traj[k] + dt[k]
    down_range[k+1] = down_range[k] + norm(x_traj[1:3,k+1] - x_traj[1:3,k])
    cross_range[k+1] = cross_range[k] + (cross(r0,v0)/norm(cross(r0,v0)))'*(x_traj[1:3,k+1] - x_traj[1:3,k])
end

bank_angle = copy(bank)
# mat"
# figure
# hold on
# plot($t_traj,$alt)
# hold off
# "
#
# mat"
# figure
# hold on
# plot($t_traj,$down_range)
# hold off
# "
#
mat"
figure
title('Bank Angle')
hold on
plot($t_traj,$bank)
hold off
"
# mat"
# figure
# title('Bank Angle Derivative')
# hold on
# stairs($t_traj(1:end-1),$σ̇)
# hold off
# "
#
mat"
figure
title('Angle of Attack')
hold on
stairs($t_traj(1:end-1),$AoA)
hold off
"
#
# mat"
# figure
# hold on
# title('dt')
# plot($dt)
# hold off
# "

r_w = 1e-2
v_w = 1e2
Q = Diagonal(SA[r_w,r_w,r_w,v_w,v_w,v_w,1e-9])
R = Diagonal(SA[1e5,1e5])

new_t, x_int, u_int, K_int, K_vec = getK(X,U,t_traj,Q,R)
@show "done"
model = EntryVehicle_fixed_time(CartesianMSLModel())
n,m = size(model)

dt = 1.0
tf = 250
t_vec = 0:dt:tf
N = length(t_vec)
X = [@SVector zeros(7) for i = 1:N]
U = [@SVector zeros(2) for i = 1:N]
tscale = 3600
using Random
Random.seed!(1);
X[1] = x0 + SVector{7}([10*normalize(randn(3));.1*tscale*normalize(randn(3));0])
# X[1] = x0
for i = 1:(N-1)

    t = t_vec[i]
    dx = X[i]  - x_int(t)
    # u = u_int(t) - K_int(t)*dx
    U[i] = u_int(t)[1:2] - K_vec[Int(floor(t))+1]*dx

    # if i == 60
    #     @infiltrate
    #     error()
    # end
    z = KnotPoint(X[i],U[i],dt)
    X[i+1] = discrete_dynamics(EntryVehicleRK_fixed_time,model,z)
end


using MATLAB
using Attitude
xm_altro = mat_from_vec(states(solver))
xm_sim = mat_from_vec(X)
mat"

figure
hold on
subplot(3,1,1)
hold on
plot($t_traj, $xm_altro(1,:))
plot($t_vec,$xm_sim(1,:))
subplot(3,1,2)
hold on
plot($t_traj, $xm_altro(2,:))
plot($t_vec,$xm_sim(2,:))
subplot(3,1,3)
hold on
plot($t_traj, $xm_altro(3,:))
plot($t_vec,$xm_sim(3,:))
hold off
"
mat"
figure
hold on
subplot(3,1,1)
hold on
plot($t_traj, $xm_altro(4,:))
plot($t_vec,$xm_sim(4,:))
subplot(3,1,2)
hold on
plot($t_traj, $xm_altro(5,:))
plot($t_vec,$xm_sim(5,:))
subplot(3,1,3)
hold on
plot($t_traj, $xm_altro(6,:))
plot($t_vec,$xm_sim(6,:))
hold off
"


um_altro = mat_from_vec(controls(solver))
um_sim = mat_from_vec(U)

mat"
figure
hold on
subplot(2,1,1)
hold
title('Bank Angle Derivative')
plot($t_traj(1:end-1),$um_altro(1,:))
plot($t_vec,$um_sim(1,:))

subplot(2,1,2)
title('Angle of Attack')
hold on
plot($t_traj(1:end-1),$um_altro(2,:))
plot($t_vec,$um_sim(2,:))

hold off
"

# X_altro = states(solver)
# rnorm = [norm(X_altro[i][1:3]) for i = 1:length(X_altro)]
# vnorm = [norm(X_altro[i][4:6]) for i = 1:length(X_altro)]
#
# mat"
# figure
# hold on
# plot($t_traj,$rnorm)
# plot($t_traj,$vnorm)
# hold off
# "
