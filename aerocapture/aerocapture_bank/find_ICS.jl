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

include(joinpath(@__DIR__,"PDIP/PDIP.jl"))
include(joinpath(@__DIR__,"dynamics.jl"))
include(joinpath(@__DIR__,"rollout_stuff.jl"))
include(joinpath(@__DIR__,"mpc.jl"))
include(joinpath(@__DIR__,"post_process.jl"))

function first_test()


alt0= 145.0
γ0 = deg2rad(-15)
V0 = 12.845*3600

# evmodel = CartesianMSLModel()
dscale = 1.0e3
tscale = 1.0
uscale = 1.0
model = EntryVehicle(CartesianMSLModel(),dscale,tscale,uscale)
#Initial conditions for MSL
Rm = model.evmodel.planet.R
r0 = [Rm+alt0, 0.0, 0.0] #Atmospheric interface at 125 km altitude
# @show norm(r0)
# error()
# V0 = 5.845*3600 #Mars-relative velocity at interface of 5.845 km/sec
# γ0 = -15.474*(pi/180.0) #Flight path angle at interface
v0 = V0*[sin(γ0), cos(γ0), 0.0]
epsilon0 = epsilon(model,[r0;v0])
σ0 = deg2rad(30)
# x0 = [r0;v0;σ0;epsilon0]
x0 = [r0/dscale;v0/(dscale/tscale); σ0; epsilon0]
# x0 = [x0;x0;x0]

# x0 = [3443.300786841311, 270.4345771068569, 0.0, -6051.64651501579, 20222.23824790719, 0.0]

# first rollout
dt = 4/3600/tscale
N = 50
X = NaN*[@SArray zeros(8) for i = 1:N]
U = [@SArray zeros(1) for i = 1:N-1]

X[1] = deepcopy(x0)
end_idx = NaN
for i = 1:(N-1)
    U[i] = [0]
    X[i+1] = rk4(model,X[i],U[i],dt)
end


# let's try some CPAG stuff
# rollout(model,x0,U_in,dt)
# -0.4 is the goal epsilon
ϵ_f = -0.4

T = 10
eps_hist = zeros(T)

althist = [zeros(N) for i = 1:T]
crhist = [zeros(N) for i = 1:T]
drhist = [zeros(N) for i = 1:T]
epshist = [zeros(N) for i = 1:T]
σhist = [zeros(N-1) for i = 1:T]

nduhist = zeros(T)
num2plot = 0
for i = 1:T

    X = rollout(model,x0,U,dt)
    althist[i],crhist[i],drhist[i] = postprocess(model,X,x0)
    epshist[i] = [X[i][8] for i = 1:N]
    σhist[i] = [X[i][7] for i = 1:N]
    eps_hist[i] = X[end][8] # 8 15 23

    # linearize
    A,B = getAB(model,X,U,dt)

    # solve cvx prob (correct)
    # U = eg_mpc(model,A,B,deepcopy(X),deepcopy(U),ϵ_f)
    U,normdu = eg_mpc_pdip(model,A,B,deepcopy(X),deepcopy(U),ϵ_f)
    nduhist[i] = normdu
    # Uhist[i] = U
    # Ut = [U[i][1] for i = 1:]
    if normdu < 1e-4
        num2plot = i
        break
    end

end
t_vec = (0:(length(X)-1))*dt*3600

# @infiltrate
# error()

# @infiltrate
# error()
pme = abs.(eps_hist .- ϵ_f)
T = num2plot
mat"
figure
hold on
%plot($eps_hist)
title('Terminal Specific Energy Error Iteration History')
plot(0:($T-1),$pme(1:$T),'linewidth',4)
set(gca, 'YScale', 'log')
xlabel('Iteration Number')
ylabel('Terminal specific energy error')
set(gca,'FontSize',16)
hold off
saveas(gcf,'eps2.eps','epsc')
"
mat"
figure
hold on
%plot($eps_hist)
title('Correction Norm Iteration History')
plot(0:($T-1),$nduhist(1:$T),'linewidth',4)
set(gca, 'YScale', 'log')
xlabel('Iteration Number')
ylabel('Terminal specific energy error')
set(gca,'FontSize',16)
hold off
saveas(gcf,'ndu.eps','epsc')
"
# mat"
# figure
# hold on
# %plot($eps_hist)
# plot(0:($T-1),$pme,'linewidth',4)
# set(gca, 'YScale', 'log')
# xlabel('Predictor-corrector iteration number')
# ylabel('Terminal specific energy error')
# set(gca,'FontSize',16)
# hold off
# saveas(gcf,'eps2.eps','epsc')
# "
# mat"
# figure
# hold on
# %plot($eps_hist)
# plot(0:($T-1),$pme,'linewidth',4)
# set(gca, 'YScale', 'log')
# xlabel('Predictor-corrector iteration number')
# ylabel('Terminal specific energy error')
# set(gca,'FontSize',16)
# hold off
# "
# eps2 = zeros(length(X))
# for i = 1:length(X)
#     eps2[i] = epsilon(model,X[i])
# end

# mat"
# figure
# hold on
# %plot($t_vec,abs($eps2 - $ϵ_f))
# plot($t_vec,$eps2)
# hold off
# "

Xm = mat_from_vec(X)
# mat"
# figure
# hold on
# %title('Bank Angle')
#
# plot($t_vec, rad2deg($Xm(7,:)),'linewidth',4)
# xlabel('Time (s)')
# ylabel('Bank Angle (degrees)')
# %set(gca,'FontSize',15)
# hold off
# "

alt = [norm(X[i][1:3])*1e3-model.evmodel.planet.R for i = 1:length(X)]

# mat"
# figure
# hold on
# title('Altitude')
# plot($t_vec, $alt,'linewidth',4)
# xlabel('Time (s)')
# ylabel('Altitude (km)')
# %set(gca,'FontSize',15)
# hold off
# "

# e_hist = [peri_apo(model,X[i]) for i = 1:length(X)]
# rp_hist = zeros(length(X))
# ra_hist = zeros(length(X))
# for i = 1:length(X)
#     rp_hist[i], ra_hist[i] = peri_apo(model,X[i])
# end
# mat"
# figure
# hold on
# plot($rp_hist)
# plot($ra_hist)
# legend('Perigee','Apogee')
# hold off
# "

# Um = mat_from_vec(U)
#
# mat"
# figure
# hold on
# plot($Um')
# hold off
# "



# --------- general plotting -----------
# alt = [norm(X[i][1:3])-model.evmodel.planet.R for i = 1:length(X)]
# #
# epsilon2 = zeros(length(X))
# for i = 1:length(X)
#     μ =  model.evmodel.planet.gravity.μ
#     r = X[i][1:3]
#     v = X[i][4:6]
#     epsilon2[i] = (dot(v,v)/2 - μ/norm(r))/1e8
# end
# # # @infiltrate
# #
# mat"
# figure
# hold on
# title('Altitude')
# plot($alt)
# hold off
# "
#
# mat"
# figure
# hold on
# title('post process epsilon')
# plot($epsilon2)
# hold off
# "
#
# Xm = mat_from_vec(X)
# eps = Xm[8,:]
# mat"
# figure
# hold on
# title('dynamics epsilon')
# plot($eps)
# hold off
# "
num2plot = float(num2plot)
# mat"
# figure
# hold on
# rgb1 = [29 38 113]/255;
# rgb2 = 1.3*[195 55 100]/255;
# drgb = rgb2-rgb1;
# cmap = [];
# for i = 1:round($num2plot)
#     px = $drhist{i};
#     py = $crhist{i};
#     plot(px,py,'Color',rgb1 + drgb*(i)/($num2plot),'linewidth',3)
#     cmap = [cmap;rgb1 + drgb*(i)/($num2plot)];
#     %plot(px(1),py(1),'r.','markersize',20)
# end
# "

# colormap(cmap);
# pp = colorbar;
# pp.Ticks = 0:(1/$num2plot):1;
# pp.Location = 'northoutside';
# pp.TickLabels = 0:round($num2plot);
# xlabel('downrange (km)')
# ylabel('crossrange (km)')
# xlim([250,700])
# %ylim([0,20])
# hold off
# %fleg = legend('figure()');
# %set(fleg,'visible','off')
# "

# this one is for plotting
mat"
figure
hold on
rgb1 = [29 38 113]/255;
rgb2 = 1.3*[195 55 100]/255;
drgb = rgb2-rgb1;
cmap = [];
for i = 1:round($num2plot)
    px = $drhist{i};
    alt = $althist{i};
    plot($t_vec,alt,'Color',rgb1 + drgb*(i)/($num2plot),'linewidth',4)
    cmap = [cmap;rgb1 + drgb*(i)/($num2plot)];
end
hold on
title('Altitude Iteration History')
colormap(cmap);
pp = colorbar;
pp.Ticks = 0:(1/$num2plot):1;
%pp.Location = 'northoutside';
pp.TickLabels = 0:round($num2plot);
%plot([0,800],ones( 2,1)*10,'r' )
%xlim([400,700])
%ylim([8,30])
xlabel('time (seconds)')
ylabel('altitude (km)')
hold off
set(gca,'FontSize',16)
saveas(gcf,'alt.eps','epsc')
"

mat"
figure
hold on
rgb1 = [29 38 113]/255;
rgb2 = 1.3*[195 55 100]/255;
drgb = rgb2-rgb1;
cmap = [];
for i = 1:round($num2plot)
    px = $drhist{i};
    alt = $epshist{i};
    plot($t_vec,alt,'Color',rgb1 + drgb*(i)/($num2plot),'linewidth',3)
    cmap = [cmap;rgb1 + drgb*(i)/($num2plot)];
end
hold on
title('Specific Orbital Energy Iteration History')
colormap(cmap);
pp = colorbar;
pp.Ticks = 0:(1/$num2plot):1;
%pp.Location = 'northoutside';
pp.TickLabels = 0:round($num2plot);
plot([0,$t_vec(end)],ones( 2,1)*$ϵ_f,'r--' ,'linewidth',2)
%xlim([400,700])
%ylim([8,30])
xlabel('time (seconds)')
ylabel('Specific Orbital Energy')
hold off
set(gca,'FontSize',16)
saveas(gcf,'eps.eps','epsc')
"

mat"
figure
hold on
rgb1 = [29 38 113]/255;
rgb2 = 1.3*[195 55 100]/255;
drgb = rgb2-rgb1;
cmap = [];
for i = 1:round($num2plot)
    px = $drhist{i};
    alt = $σhist{i};
    plot($t_vec,rad2deg(alt),'Color',rgb1 + drgb*(i)/($num2plot),'linewidth',3)
    cmap = [cmap;rgb1 + drgb*(i)/($num2plot)];
end
hold on
title('Bank Angle Iteration History')
colormap(cmap);
pp = colorbar;
pp.Ticks = 0:(1/$num2plot):1;
%pp.Location = 'northoutside';
pp.TickLabels = 0:round($num2plot);
%plot([0,$t_vec(end)],ones( 2,1)*$ϵ_f,'r--' ,'linewidth',2)
%xlim([400,700])
%ylim([8,30])
xlabel('time (seconds)')
ylabel('Bank Angle (deg)')
hold off
set(gca,'FontSize',16)
saveas(gcf,'bank.eps','epsc')
"
    return nothing
end
first_test()



# mat"
# figure
# hold on
# plot($xm(1:3,:)')
# hold off
# "
#
# mat"
# figure
# hold on
# plot($tt,$xm(3,:),'o')
# plot($tt,$z2)
# hold off
# "


# aa = [randn(3,10) for i = 1:3]
#
# mat"
# figure
# hold on
# for i = 1:length($aa)
#     p = $aa{1}
#     disp(p)
# end
# hold off
# "
#
