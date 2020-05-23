using LinearAlgebra
using DifferentialEquations
using EntryGuidance
using Test

function open_loop_dynamics!(dx,x,model,t)
    a = [10.0*pi/180.0; -45.0*pi/180] # [AoA; Bank Angle]
    u = angles_input(a,x,model)
    dynamics!(dx,x,u,model)
end

#Initial conditions for a low-Earth orbit
planet = Earth()
vehicle = SimpleMSLVehicle()
m_v = VinhModel(planet, vehicle)
m_c = CartesianModel(planet, vehicle)
μ = planet.gravity.μ
Re = planet.R
r0 = [Re+200.0, 0, 0] #+ 10*randn(3)
v0 = [0, sqrt(μ/norm(r0)), 0] #+ 10*randn(3)
x0_cart = inertial_to_planet_fixed([r0; v0], planet)
tspan = (0.0, 1.5)

prob1 = ODEProblem(open_loop_dynamics!,x0_cart,tspan,m_c)
soln1 = solve(prob1, Tsit5(), reltol=1e-9, abstol=1e-9)

x0_vinh = cartesian_to_vinh(x0_cart)

prob2 = ODEProblem(open_loop_dynamics!,x0_vinh,tspan,m_v)
soln2 = solve(prob2, Tsit5(), reltol=1e-12, abstol=1e-12)


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
    x_2[:,k] = soln2(tsteps[k])
    x_3[:,k] = vinh_to_cartesian(soln2(tsteps[k]))
    x_4[:,k] = cartesian_to_vinh(soln1(tsteps[k]))

end

@test x_1≈x_3 atol=1e-4
@test x_2[[1,4],:]≈x_4[[1,4],:] atol=1e-4
@test mod.(x_2[[2,3,5,6],:],2π)≈mod.(x_4[[2,3,5,6],:],2π) atol=1e-4
