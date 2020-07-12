using LinearAlgebra
using StaticArrays
using EntryGuidance
const EG = EntryGuidance
using TrajectoryOptimization
const TO = TrajectoryOptimization
using RobotDynamics
const RD = RobotDynamics
using Altro

struct EntryVehicleModel{T} <: TO.AbstractModel
    evmodel::EG.CartesianModel{T}
end

function RD.dynamics(model::EntryVehicleModel, x, u)
    α̇ = u[1]-u[2]
    σ̇ = u[3]-u[4]
    [EG.dynamics(x[1:6], EG.angles_input(x[7:8],x[1:6],model.evmodel), model.evmodel); α̇; σ̇]
end

Base.size(::EntryVehicleModel) = 8,4

model = EntryVehicleModel(CartesianMSLModel())
n,m = size(model)
tf = 6/60.0
N = 121

#Initial conditions for MSL
Rm = model.evmodel.planet.R
r0 = [Rm+125.0, 0.0, 0.0] #Atmospheric interface at 125 km altitude
V0 = 5.845*3600 #Mars-relative velocity at interface of 5.845 km/sec
γ0 = -15.474*(pi/180.0) #Flight path angle at interface
v0 = V0*[sin(γ0), cos(γ0), 0.0]
α = 15.0*pi/180 #derived from MSL initial L/D = 0.24
σ = 0.0*pi/180 #this is totally made up
x0 = SVector{n}([r0; v0; α; σ])

#Final conditions for MSL 631.979 km down range and 7.869 km cross-range from entry point
rf = Rm.*[cos(7.869/Rm)*cos(631.979/Rm); cos(7.869/Rm)*sin(631.979/Rm); sin(7.869/Rm)]
vf = zeros(3)
xf = SVector{n}([rf; vf; α; σ])

#Cost function
Q = Diagonal(SVector{n}(1e-4.*ones(n)))
R = Diagonal(SVector{m}(1e-4.*ones(m)))
H = zeros(m,n)
q = -Q*xf
r = 10.0*ones(m)
c = 0.0
stage_cost = QuadraticCost(Q,R,H,q,r,c,terminal=false)

Qn = Diagonal(@SVector [100.0, 100.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0])
terminal_cost = LQRCost(Qn,R,xf,terminal=true)

obj = Objective(stage_cost,terminal_cost,N)

#Constraints
cons = TO.ConstraintList(n,m,N)
∞ = Inf64
add_constraint!(cons, BoundConstraint(n,m,x_min=SA[-∞,-∞,-∞,-∞,-∞,-∞,0.0,-pi],x_max=SA[∞,∞,∞,∞,∞,∞,20*pi/180,pi]),1:N)
add_constraint!(cons, BoundConstraint(n,m,u_min=[0.0,0.0,0.0,0.0],u_max=[∞,∞,∞,∞]),1:(N-1))

prob = TO.Problem(model, obj, xf, tf, x0=x0, constraints=cons)
solver = AugmentedLagrangianSolver(prob)
solve!(solver)

X = states(solver)
U = controls(solver)

max_violation(solver)

x_traj = zeros(n,N)
alt = zeros(N)
AoA = zeros(N)
bank = zeros(N)
for k = 1:N
    x_traj[:,k] .= X[k]
    alt[k] = norm(x_traj[1:3,k])-Rm
    AoA[k] = x_traj[7,k]
    bank[k] = x_traj[8,k]
end

u_traj = zeros(m,N-1)
α̇ = zeros(N-1)
σ̇ = zeros(N-1)
for k = 1:(N-1)
    u_traj[:,k] .= U[k]
    α̇[k] = u_traj[1,k]-u_traj[2,k]
    σ̇[k] = u_traj[3,k]-u_traj[4,k]
end

using Plots
plotlyjs()
plot(alt)
plot(AoA)
plot(bank)
plot(α̇)
plot(σ̇)
