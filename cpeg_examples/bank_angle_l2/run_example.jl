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
using Random
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
Random.seed!(123)

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
σ0 = deg2rad(3)
x0 = [r0;v0;σ0]

#Final conditions for MSL 631.979 km down range and 7.869 km cross-range from entry point
Rf = Rm+10.0 #Parachute opens at 10 km altitude
rf = Rf*[cos(7.869/Rf)*cos(631.979/Rf); cos(7.869/Rf)*sin(631.979/Rf); sin(7.869/Rf)]
vf = zeros(3)
xf = [rf;vf;0]

# scale units
x0 = [r0/dscale;v0/(dscale/tscale); σ0]
xf = [rf/dscale;vf/(dscale/tscale); σ0]

# initial control
dt = (5/3600)/tscale
N = 180
X = NaN*[@SArray zeros(7) for i = 1:N]
U = [[0.00*randn()]  for i = 1:N-1]

# number of iterations
T = 15

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
    Xc, U, dunorm[i] = eg_mpc2(model,A,B,X,U,xf)
    # Xc, U, dunorm[i] = eg_mpc_l1(model,A,B,X,U,xf)

end

# AoA, bank = processU(model::EntryVehicle,X,U)
N = length(X)
T_vec = 0:(dt*3600):((N-1)*(dt*3600))
xf_dr, xf_cr = rangedistances(model,xf,x0)

bank = [X[i][7] for i = 1:(length(X)-1)]


# # number of trajectories to plot (this has to be a float for some reason)
num2plot = float(T)
mat"
figure
hold on
%title('Bank Angle')
plot($T_vec(1:end-1)/60, rad2deg($bank),'linewidth',3)
ylabel('Bank Angle (deg)')
xlabel('time (min)')
hold off
"
#
mat"
figure
hold on
title('Du norm')
plot($dunorm)
hold off
"


plot_groundtracks(drhist,crhist,althist,xf_dr,xf_cr,num2plot)
#
# this one is for plotting
# mat"
# figure
# hold on
# rgb1 = [29 38 113]/255;
# rgb2 = 1.3*[195 55 100]/255;
# drgb = rgb2-rgb1;
# for i = 1:length($drhist)
#     px = $drhist{i};
#     py = $crhist{i};
#     if i < ($num2plot +1)
#         plot(px,py,'Color',rgb1 + drgb*(i-1)/($num2plot),'linewidth',3)
#     end
#     plot(px(1),py(1),'r.','markersize',20)
# end
# plot($xf_dr,$xf_cr,'g.','markersize',20)
# xlabel('downrange (km)')
# ylabel('crossrange (km)')
# hold off
# %saveas(gcf,'range.png')
# addpath('/Users/kevintracy/devel/WiggleSat/matlab2tikz-master/src')
# %matlab2tikz('bankaoa_track.tex')
# %close all
# "
#
# # this one is for plotting
# mat"
# figure
# hold on
# rgb1 = [29 38 113]/255;
# rgb2 = 1.3*[195 55 100]/255;
# drgb = rgb2-rgb1;
# for i = 1:length($althist)
#     px = $drhist{i};
#     alt = $althist{i};
#     if i < ($num2plot +1)
#         colo = drgb*(i-1)/$num2plot;
#         plot(px,alt,'Color',rgb1 + colo,'linewidth',3)
#     end
#     plot(px(1),alt(1),'r.','markersize',20)
# end
# plot([0,800],ones( 2,1)*10,'r' )
# plot($xf_dr,10,'g.','markersize',20)
# xlim([0 650])
# xlabel('downrange (km)')
# ylabel('altitude (km)')
# hold off
# %saveas(gcf,'alt.png')
# addpath('/Users/kevintracy/devel/WiggleSat/matlab2tikz-master/src')
# %matlab2tikz('bankaoa_alt.tex')
# %close all
# "

    return 0
end
Xsim2 = first_test()
