using LinearAlgebra
using DifferentialEquations
using EntryGuidance
using Test

function closed_loop_dynamics1!(dx,x,params,t)
    vehicle, planet = params
    u = cartesian_test_controller(x,vehicle,planet)
    cartesian_dynamics!(dx,x,u,planet)
end

function closed_loop_dynamics2!(dx,x,params,t)
    vehicle, planet = params
    u = vinh_test_controller(x,vehicle,planet)
    vinh_dynamics!(dx,x,u,planet)
end

function closed_loop_dynamics3!(dx,x,params,t)
    vehicle, planet, y = params
    u = linear_test_controller(x,vehicle,planet)
    linear_dynamics!(dx,x,u,y,planet)
end

#Initial conditions for a low-Earth orbit
planet = EarthPlanetModel()
vehicle = SimpleSphereConeVehicle()
μ = planet.gravity.μ
Re = planet.R
r0 = [Re+200.0, 0, 0] #+ 10*randn(3)
v0 = [0, sqrt(μ/norm(r0)), 0] #+ 10*randn(3)
x0_cart = inertial_to_planet_fixed([r0; v0], planet)
tspan = (0.0, 1.5)

kg = norm(gravitational_acceleration(r0,planet))/norm(r0)
y = [kg]
params = (vehicle, planet, y)

prob1 = ODEProblem(closed_loop_dynamics1!,x0_cart,tspan,params)
soln1 = solve(prob1, Tsit5(), reltol=1e-9, abstol=1e-9)

x0_vinh = cartesian_to_vinh(x0_cart)

prob2 = ODEProblem(closed_loop_dynamics2!,x0_vinh,tspan,params)
soln2 = solve(prob2, Tsit5(), reltol=1e-12, abstol=1e-12)

prob3 = ODEProblem(closed_loop_dynamics3!,x0_cart,tspan,params)
soln3 = solve(prob3, Tsit5(), reltol=1e-9, abstol=1e-9)

using Plots
plot(soln1,vars=(1,2))
plot(soln3,vars=(1,2))

#Take samples along trajectory
nsteps = 20
tsteps = LinRange(tspan[1],tspan[2],nsteps)
x_1 = zeros(6,nsteps)
x_2 = zeros(6,nsteps)
x_3 = zeros(6,nsteps)
x_4 = zeros(6,nsteps)
x_5 = zeros(6,nsteps)
for k = 1:nsteps
    x_1[:,k] = soln1(tsteps[k])
    x_2[:,k] = vinh_to_cartesian(soln2(tsteps[k]))
    x_3[:,k] = cartesian_to_vinh(soln1(tsteps[k]))
    x_4[:,k] = soln2(tsteps[k])
    x_5[:,k] = soln3(tsteps[k])
end



@test x_1≈x_2 atol=1e-4
@test x_3[[1,4],:]≈x_4[[1,4],:] atol=1e-4
@test mod.(x_3[[2,3,5,6],:],2π)≈mod.(x_4[[2,3,5,6],:],2π) atol=1e-4
