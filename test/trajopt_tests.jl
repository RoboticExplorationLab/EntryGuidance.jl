using LinearAlgebra
using StaticArrays
using EntryGuidance
const EG = EntryGuidance
using TrajectoryOptimization
const TO = TrajectoryOptimization

struct EntryVehicleModel{T} <: TO.AbstractModel
    evmodel::EG.CartesianModel{T}
end

Base.size(::EntryVehicleModel) = 8,2

function TO.dynamics(model::EntryVehicleModel, x, u)
    [EG.dynamics(x[1:6], EG.angles_input(x[7:8],x[1:6],model.evmodel), model.evmodel); u]
end

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

Q = Diagonal(SVector{n}(1e-6.*ones(n)))
Qn = Diagonal(@SVector [100.0, 100.0, 100.0, 0.1, 0.1, 0.1, 0.0, 0.0])
R = Diagonal(SVector{m}(0.1.*ones(m)))
obj = LQRObjective(Q,R,Qn,xf,N)

cons = TO.ConstraintSet(n,m,N)
∞ = Inf64
add_constraint!(cons, BoundConstraint(n,m,x_min=SA[-∞,-∞,-∞,-∞,-∞,-∞,0.0,-pi],x_max=SA[∞,∞,∞,∞,∞,∞,20*pi/180,pi]), 1:N)

prob = TO.Problem(model, obj, xf, tf, x0=x0, constraints=cons)
solver = AugmentedLagrangianSolver(prob)
solve!(solver)

X = states(solver)
U = controls(solver)

using Plots

plot(X)
