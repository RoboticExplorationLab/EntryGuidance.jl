using LinearAlgebra
using DifferentialEquations
using EntryGuidance
using Test

function open_loop_dynamics!(dx,x,params,t)
    model, α, σ = params
    u = angles_input([α, σ],x,model)
    dynamics!(dx,x,u,model)
end

#Rough initial conditions for MSL from
#https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20130010129.pdf
model = CartesianMSLModel()
Rm = model.planet.R
r0 = [Rm+125.0, 0.0, 0.0] #Atmospheric interface at 125 km altitude, 631.979 km down range and 7.869 km cross-range of target
V0 = 5.845*3600 #Mars-relative velocity at interface of 5.845 km/sec
γ0 = -15.474*(pi/180.0) #Flight path angle at interface
v0 = V0*[sin(γ0), cos(γ0), 0.0]
x0 = [r0; v0]
α = 15.0*pi/180 #derived from MSL initial L/D = 0.24
σ = 0.0*pi/180 #this is totally made up
params = (model, α, σ)

tspan = (0.0, 0.1)

prob = ODEProblem(open_loop_dynamics!,x0,tspan,params)
soln = solve(prob, Tsit5(), reltol=1e-9, abstol=1e-9)

using Plots
# plotly()
# gr()

nsteps = 100
t = LinRange(soln.t[1],soln.t[end],nsteps)
alt = zeros(nsteps)
for k = 1:nsteps
    alt[k] = norm(soln(t[k])[1:3])-Rm
end
plot(t.*60,alt)
