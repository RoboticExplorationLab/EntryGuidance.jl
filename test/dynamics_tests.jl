using LinearAlgebra
#using Test
using DifferentialEquations
using EntryGuidance

function closed_loop_dynamics1!(dx,x,params,t)
    vehicle, planet = params
    u = cartesian_no_lift(x,vehicle,planet)
    cartesian_dynamics!(dx,x,u,planet)
end

function closed_loop_dynamics2!(dx,x,params,t)
    vehicle, planet = params
    u = vinh_no_lift(x,vehicle,planet)
    vinh_dynamics!(dx,x,u,planet)
end

#Initial conditions for a low-Earth orbit
planet = EarthPlanetModel()
vehicle = SimpleSphereConeVehicle()
params = (vehicle, planet)
μ = planet.gravity.μ
Re = planet.R
r0 = [Re+200.0, 0, 0] + 10*randn(3)
v0 = [0, sqrt(μ/norm(r0)), 0] + 10*randn(3)
x0_cart = inertial_to_planet_fixed([r0; v0], planet)
tspan = (0.0, 1.5)

prob1 = ODEProblem(closed_loop_dynamics1!,x0_cart,tspan,params)
soln1 = solve(prob1, Tsit5(), reltol=1e-9, abstol=1e-9)

x0_vinh = cartesian_to_vinh(x0_cart)

prob2 = ODEProblem(closed_loop_dynamics2!,x0_vinh,tspan,params)
soln2 = solve(prob2, Tsit5(), reltol=1e-14, abstol=1e-14)

tsteps = LinRange(tspan[1],tspan[2],100)
x_1 = zeros(6,100)
x_2 = zeros(6,100)
x_3 = zeros(6,100)
x_4 = zeros(6,100)
for k = 1:100
    x_2[:,k] = vinh_to_cartesian(soln2(tsteps[k]))
    x_1[:,k] = soln1(tsteps[k])
    x_3[:,k] = cartesian_to_vinh(soln1(tsteps[k]))
    x_4[:,k] = soln2(tsteps[k])
end

using Plots
plotlyjs()
plot(soln,vars=(1,2))
