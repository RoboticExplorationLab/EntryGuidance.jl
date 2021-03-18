using RobotDynamics
const RD = RobotDynamics
using TrajectoryOptimization
const TO = TrajectoryOptimization
using Attitude
using MATLAB
using Infiltrator
using StaticArrays
using JLD2
using ForwardDiff
using EntryGuidance
const EG = EntryGuidance
@load "goodtraj.jld2" traj
struct EntryVehicle{T} <: TO.AbstractModel
    evmodel::EG.CartesianModel{T}
end
function evdynamics(model::EntryVehicle, x, u)
    #unpack control
    α = u[1] #angle of attack
    σ = u[2] #bank angle

    #unpack state
    r = x[1:3]
    v = x[4:6]
    V = norm(v)

    #atmospheric density
    ρ = atmospheric_density(r, model.evmodel)

    #Calculate drag acceleration
    Cd = drag_coefficient(α, model.evmodel)
    A = model.evmodel.vehicle.A
    m = model.evmodel.vehicle.m
    D = 0.5*Cd*ρ*A*V*V/m

    #Calculate lift acceleration
    Cl = lift_coefficient(α, model.evmodel)
    L = 0.5*Cl*ρ*A*V*V/m

    #get gravity
    g = gravitational_acceleration(r, model.evmodel)

    #Terms involving Ω (planet rotation) are often thrown out in the literature.
    Ω = model.evmodel.planet.Ω #Set Ω = 0.0 here if you want that behavior
    Ω̂ = hat([0, 0, Ω])

    #Aerodynamic acceleration
    e1 = normalize(cross(r,v))
    e2 = normalize(cross(v,e1))
    a = -(D/norm(v))*v + L*sin(σ)*e1 + L*cos(σ)*e2

    v̇ = a + g - 2*Ω̂*v - Ω̂*Ω̂*r

    return [v; v̇]
end
function lin_evdynamics(model::EntryVehicle, x, u)

    α = 10.0

    #unpack state
    r = x[1:3]
    v = x[4:6]
    V = norm(v)

    #atmospheric density
    ρ = atmospheric_density(r, model.evmodel)

    #Calculate drag acceleration
    Cd = drag_coefficient(α, model.evmodel)
    A = model.evmodel.vehicle.A
    m = model.evmodel.vehicle.m
    D = 0.5*Cd*ρ*A*V*V/m

    #Calculate lift acceleration
    # Cl = lift_coefficient(α, model.evmodel)
    # L = 0.5*Cl*ρ*A*V*V/m

    #get gravity
    g = gravitational_acceleration(r, model.evmodel)

    #Terms involving Ω (planet rotation) are often thrown out in the literature.
    Ω = model.evmodel.planet.Ω #Set Ω = 0.0 here if you want that behavior
    Ω̂ = hat([0, 0, Ω])

    #Aerodynamic acceleration
    e1 = normalize(cross(r,v))
    e2 = normalize(cross(v,e1))
    # a = -(D/norm(v))*v + L*sin(σ)*e1 + L*cos(σ)*e2
    a = -(D/norm(v))*v + u[1]*e1 + u[2]*e2

    # replace a with u

    v̇ = a + g - 2*Ω̂*v - Ω̂*Ω̂*r

    return [v; v̇]
end
# Define the discrete dynamics function
function rk3(model::EntryVehicle,
        x::StaticVector, u::StaticVector, t, dt)

    h = dt/3600.0 #u is in seconds, dynamics are in hours

    k1 = evdynamics(model, x,             u)*h;
    k2 = evdynamics(model, x + k1/2,      u)*h;
    k3 = evdynamics(model, x - k1 + 2*k2, u)*h;

    return x + (k1 + 4*k2 + k3)/6
end
Base.size(::EntryVehicle) = 6,2
model = EntryVehicle(CartesianMSLModel())
n,m = size(model)

dt = 0.1
N = length(traj.t)

X = [@SVector zeros(n) for i = 1:N]
U = [@SVector zeros(m) for i = 1:N-1]
X[1] = traj.X[1][1:6]

for i = 1:N-1

    α = traj.U[i][2]
    α = deg2rad(10)
    σ = traj.X[i][7]
    U[i] = SA[α,σ]
    t = traj.t[i]
    X[i+1] = rk3(model,X[i],U[i],t,dt)

end


Um = mat_from_vec(U)

mat"
figure
hold on
plot(rad2deg($Um'))
legend('Angle of Attack','Bank Angle')
ylabel('Degrees')
hold off
"

alt = [(norm(X[i][1:3]) - model.evmodel.planet.R) for i = 1:length(X)]

mat"
figure
hold on
plot($alt)
hold off
"

xm_sim = mat_from_vec(X)
t_vec = traj.t
mat"

figure
hold on
sgtitle('MCMF Position')
subplot(3,1,1)
hold on
plot($t_vec,$xm_sim(1,:))
subplot(3,1,2)
hold on
plot($t_vec,$xm_sim(2,:))
subplot(3,1,3)
hold on
plot($t_vec,$xm_sim(3,:))
hold off
"
mat"
figure
sgtitle('MCMF Velocity')
hold on
subplot(3,1,1)
hold on
plot($t_vec,$xm_sim(4,:))
subplot(3,1,2)
hold on
plot($t_vec,$xm_sim(5,:))
subplot(3,1,3)
hold on
plot($t_vec,$xm_sim(6,:))
hold off
"
