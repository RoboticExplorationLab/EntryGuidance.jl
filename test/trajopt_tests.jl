using LinearAlgebra
using StaticArrays
using EntryGuidance
const EG = EntryGuidance
using TrajectoryOptimization
const TO = TrajectoryOptimization
using RobotDynamics
const RD = RobotDynamics
using Altro

struct EntryVehicle{T} <: TO.AbstractModel
    evmodel::EG.CartesianModel{T}
end

function RD.dynamics(model::EntryVehicle, x, u)
    α = 15.0*pi/180
    σ̇ = u[1]-u[2]
    [EG.dynamics(x[1:6], EG.angles_input([α, x[7]],x[1:6],model.evmodel), model.evmodel); σ̇; 0.0]
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

    return [x[1:7] + (k1[1:7] + 4*k2[1:7] + k3[1:7])/6; x[8]+u[3]]
end

Base.size(::EntryVehicle) = 8,3

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
α = 15.0*pi/180 #derived from MSL initial L/D = 0.24
σ = 0.0*pi/180 #this is totally made up
x0 = SVector{n}([r0; v0; σ; 0.0])

#Final conditions for MSL 631.979 km down range and 7.869 km cross-range from entry point
Rf = Rm+10.0 #Parachute opens at 10 km altitude
rf = Rf*[cos(7.869/Rf)*cos(631.979/Rf); cos(7.869/Rf)*sin(631.979/Rf); sin(7.869/Rf)]
vf = zeros(3)
xf = SVector{n}([rf; vf; σ; tf*3600])

#Cost function
Q = Diagonal(SVector{n}(zeros(n)))
R = Diagonal(SVector{m}(1e-5.*ones(m)))
H = zeros(m,n)
q = -Q*xf
r = SA[10.0, 10.0, 0.0]
c = 0.0
stage_cost = DiagonalCost(Q,R,H,q,r,c,terminal=false)

Qn = Diagonal(@SVector [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1e-5])
terminal_cost = LQRCost(Qn,R,xf,terminal=true)

obj = Objective(stage_cost,terminal_cost,N)

#Constraints
cons = TO.ConstraintList(n,m,N)
∞ = Inf64
add_constraint!(cons, BoundConstraint(n,m,x_min=SA[-∞,-∞,-∞,-∞,-∞,-∞,-pi,-∞],x_max=SA[∞,∞,∞,∞,∞,∞,pi,∞]),1:N)
add_constraint!(cons, BoundConstraint(n,m,u_min=[0.0,0.0,1.0],u_max=[∞,∞,4.0]),1:(N-1))
add_constraint!(cons, GoalConstraint(xf, [1,2,3]), N:N)

#Initial Controls
U0 = abs.(randn(m,N-1))
U0[3,:] .= dt0*ones(N-1)

prob = TO.Problem(model, obj, xf, tf, x0=x0, U0=U0, constraints=cons, integration=EntryVehicleRK)
solver = AugmentedLagrangianSolver(prob)
solve!(solver)

X = states(solver)
U = controls(solver)

max_violation(solver)

#Test rollout
x_traj = zeros(n,N)
# x_traj[:,1] .= x0
# for k = 1:(N-1)
#     x_traj[:,k+1] .= discrete_dynamics(EntryVehicleRK,model,SVector{n}(x_traj[:,k]),SVector{m}(U0[:,k]),0.0,2.0)
# end
alt = zeros(N)
#AoA = zeros(N)
bank = zeros(N)
ts = zeros(N)
for k = 1:N
    x_traj[:,k] .= X[k]
    alt[k] = norm(x_traj[1:3,k])-Rm
    #AoA[k] = x_traj[7,k]
    bank[k] = x_traj[7,k]
    ts[k] = x_traj[8,k]
end

u_traj = zeros(m,N-1)
# u_traj = U0
σ̇ = zeros(N-1)
dt = zeros(N-1)
for k = 1:(N-1)
    u_traj[:,k] .= U[k]
    σ̇[k] = u_traj[1,k]-u_traj[2,k]
    dt[k] = u_traj[3,k]
end

using Plots
#plotlyjs()
pyplot()
plot(alt)
#plot(AoA)
plot(bank)
#plot(α̇)
plot(σ̇)
plot(dt)
plot(ts)
