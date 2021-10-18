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
using JLD2
using SuiteSparse
using SparseArrays



include(joinpath(@__DIR__,"dynamics.jl"))
include(joinpath(@__DIR__,"rollout_stuff.jl"))
include(joinpath(@__DIR__,"mpc.jl"))
include(joinpath(@__DIR__,"post_process.jl"))
include(joinpath(dirname(@__DIR__),"plotting_recipes.jl"))

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
dt = (2/3600)/tscale
N = 180
X = NaN*[@SArray zeros(6) for i = 1:N]
U = [[-0.05;0.55] + 0.00*randn(2) for i = 1:N-1]

# number of iterations
T = 7

# vectors for storing trajectory information
althist = [zeros(2) for i = 1:T]
drhist = [zeros(2) for i = 1:T]
crhist = [zeros(2) for i = 1:T]
dunorm = zeros(T)

# main loop
for i = 1:T

    # prediction
    X, U = rollout(model,deepcopy(x0),U,dt)

    # pull out altitude, downrange, and crossrange paths for plotting
    althist[i], drhist[i], crhist[i] = postprocess(model,X,x0)

    # jacobians (linearization)
    A,B = getAB(model,X,U,dt)

    # correction (convex solve)
    Xc, U, dunorm[i] = eg_mpc(model,A,B,X,U,xf,i)

end


AoA, bank = processU(model::EntryVehicle,X,U)
N = length(X)
T_vec = 0:(dt*3600):((N-1)*(dt*3600))
xf_dr, xf_cr = rangedistances(model,xf,x0)


# number of trajectories to plot (this has to be a float for some reason)
num2plot = float(T)
fl_traj = (T_vec = T_vec, bank = bank, AoA = AoA,
           alt = althist[end], dr = drhist[end],cr = crhist[end])

# jldsave("cpeg_examples/bank_angle/trajectories/full_lift.jld2";fl_traj)
flb = fl_traj.bank
fla = fl_traj.AoA
flt = fl_traj.T_vec

mat"
figure
hold on
plot($flt(1:end-1)/60, rad2deg($flb),'linewidth',3)
plot($flt(1:end-1)/60, rad2deg($fla),'linewidth',3)
legend('bank angle','angle of attack','location','northwest')
ylabel('aero angle (deg)')
xlabel('time (min)')
%addpath('/Users/kevintracy/devel/WiggleSat/matlab2tikz-master/src')
%matlab2tikz(strcat('cpeg_examples/bank_angle/tikz/fl_controls.tex'))
"

# number of trajectories to plot (this has to be a float for some reason)
num2plot = 4.0#float(T)
plot_groundtracks(drhist,crhist,althist,xf_dr,xf_cr,num2plot,"FL")

    return 0
end
Xsim2 = first_test()
