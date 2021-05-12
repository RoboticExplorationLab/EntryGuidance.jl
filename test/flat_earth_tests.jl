using EntryGuidance
using LinearAlgebra
using StaticArrays
using MATLAB
using Attitude

struct EntryVehicle{T} <: TO.AbstractModel
    evmodel::EG.CartesianModel{T}
end


function evdynamics(model::EntryVehicle, x, u,atmo_uncert)
    #unpack control
    # σ̇ = u[1]
    # α = u[2] #angle of attack
    # σ = x[7] #bank angle
    α = deg2rad(18)
    #unpack state
    r = x[1:3]
    v = x[4:6]
    V = norm(v)

    #atmospheric density
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
function rk4(model,x_n,u,dt,atmo_uncert)
    k1 = dt*evdynamics(model,x_n,u,atmo_uncert)
    k2 = dt*evdynamics(model,x_n+k1/2,u,atmo_uncert)
    k3 = dt*evdynamics(model,x_n+k2/2,u,atmo_uncert)
    k4 = dt*evdynamics(model,x_n+k3,u,atmo_uncert)
    return (x_n + (1/6)*(k1 + 2*k2 + 2*k3 + k4))
end
struct EntryVehicle{T} <: TO.AbstractModel
    evmodel::EG.CartesianModel{T}
end
model = EntryVehicle(CartesianMSLModel())
# n,m = size(model)
n = 6
m = 3
tf = 5.5/60.0 #time is in hours
N = 121
dt0 = (tf*3600)/(N-1)

#Initial conditions for MSL
Rm = model.evmodel.planet.R
r0 = [125.0, 0.0, 0.0] #Atmospheric interface at 125 km altitude
V0 = 5.845*3600 #Mars-relative velocity at interface of 5.845 km/sec
γ0 = -15.474*(pi/180.0) #Flight path angle at interface
v0 = V0*[sin(γ0), cos(γ0), 0.0]
α0 = 15.0*pi/180 #derived from MSL initial L/D = 0.24
σ = 0.0*pi/180 #this is totally made up
x0 = SVector{n}([r0; v0])

N = 80
X = [zeros(6) for i = 1:N]
U = [[0;0.8] for i =1:N-1]
X[1] = deepcopy(x0)

for i = 1:(N-1)
    X[i+1] = rk4(model,X[i],U[i],2/3600,1)
end

xm = mat_from_vec(X)
mat"
figure
hold on
plot($xm(2,:),$xm(3,:))
hold off
"

mat"
figure
hold on
plot($xm(1,:))
hold off
"
# mat"
# figure
# hold on
# plot($xm(7,:))
# hold off
# "
um = mat_from_vec(U)

mat"
figure
hold on
plot($um')
hold off
"
