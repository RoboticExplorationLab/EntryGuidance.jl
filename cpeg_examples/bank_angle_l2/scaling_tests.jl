using LinearAlgebra
using EntryGuidance
const EG = EntryGuidance
using StaticArrays
using Attitude
using MATLAB
using Infiltrator
using ForwardDiff
using Convex
using Mosek, MosekTools
# using COSMO
using SuiteSparse
using SparseArrays
# using Interpolations


include(joinpath(@__DIR__,"dynamics.jl"))
include(joinpath(@__DIR__,"rollout_stuff.jl"))
include(joinpath(@__DIR__,"mpc.jl"))
include(joinpath(@__DIR__,"post_process.jl"))


function dist_from_target(X,xf)
    rr = normalize(xf[1:3])
    Qn = I - rr*rr'
    return norm(Qn*(X[end][1:3]-xf[1:3]))
end

function first_test()
dscale = 1.0e3
tscale = 1.0
uscale = 1.0
model = EntryVehicle(CartesianMSLModel(),dscale,tscale,uscale)
# @infiltrate
# error()
#Initial conditions for MSL
Rm = model.evmodel.planet.R
r0 = [Rm+125.0, 0.0, 0.0] #Atmospheric interface at 125 km altitude
V0 = 5.845*3600 #Mars-relative velocity at interface of 5.845 km/sec
γ0 = -15.474*(pi/180.0) #Flight path angle at interface
v0 = V0*[sin(γ0), cos(γ0), 0.0]
x0 = [r0;v0]

#Final conditions for MSL 631.979 km down range and 7.869 km cross-range from entry point
Rf = Rm+10.0 #Parachute opens at 10 km altitude
rf = Rf*[cos(7.869/Rf)*cos(631.979/Rf); cos(7.869/Rf)*sin(631.979/Rf); sin(7.869/Rf)]
vf = zeros(3)
xf = [rf;vf]

# scale units
x0 = [r0/dscale;v0/(dscale/tscale)]
xf = [rf/dscale;vf/(dscale/tscale)]

# first rollout to 10km altitude
dt = (2/3600)/tscale # scaled
N = 180
X = NaN*[@SArray zeros(6) for i = 1:N]
U = [[0.03;0.5] for i = 1:N-1]

T = 7
althist = [zeros(2) for i = 1:T]
drhist = [zeros(2) for i = 1:T]
crhist = [zeros(2) for i = 1:T]
dunorm = zeros(T)


X, U, t_vec, t_impact = rollout(model,deepcopy(x0),U,dt)


rnorm = [norm(X[i][1:3]) for i = 1:length(X)]
vnorm = [norm(X[i][4:6]) for i = 1:length(X)]
#
# Xs = deepcopy(X);
# dscale = 1e3
# tscale = 1
# for i = 1:length(X)
#     r = X[i][1:3]
#     v = X[i][4:6]
#
#     rs = r/dscale
#     vs = v/(dscale/tscale)
#
#     Xs[i] = [rs;vs]
# end
#
# rsnorm = [norm(Xs[i][1:3]) for i = 1:length(Xs)]
# vsnorm = [norm(Xs[i][4:6]) for i = 1:length(Xs)]
#
#
# mat"
# figure
# hold on
# title('scaled')
# plot($rsnorm)
# plot($vsnorm)
# "

mat"
figure
hold on
title('unscaled')
plot($rnorm)
plot($vnorm)
"

return nothing
end

first_test()
