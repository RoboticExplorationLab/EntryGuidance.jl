using LinearAlgebra
using Plots
using Test
using DifferentialEquations
using EntryGuidance

function closed_loop_dynamics!(dx,x,params,t)
    vehicle, planet = params
    u = cartesian_no_lift(x,vehicle,planet)
    cartesian_dynamics!(dx,x,u,planet)
end

#Initial conditions for a low-Earth orbit
planet = EarthPlanetModel()
vehicle = SimpleSphereConeVehicle()
params = (vehicle, planet)
μ = planet.gravity.μ
Re = planet.R
r0 = [Re+270.0, 0, 0] #+ 10*randn(3)
v0 = [0, sqrt(μ/norm(r0)), 0] #+ 0.01*randn(3)
x0 = inertial_to_planet_fixed([r0; v0], planet)
tspan = (0.0, 10*90*60.0)

prob = ODEProblem(closed_loop_dynamics!,x0,tspan,params)
soln = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)

plot(soln,vars=(1,2))
