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
x0 = [x0;x0;x0]

# x0 = [3443.300786841311, 270.4345771068569, 0.0, -6051.64651501579, 20222.23824790719, 0.0]

# first rollout
dt = 4/3600/tscale
N = 70
X = NaN*[@SArray zeros(8*3) for i = 1:N]
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

T = 20
eps_hist = [zeros(3) for i = 1:T]
for i = 1:T

    X = rollout(model,x0,U,dt)
    eps_hist[i] = [X[end][8], X[end][16], X[end][24]] # 8 15 23

    # linearize
    A,B = getAB(model,X,U,dt)

    # solve cvx prob (correct)
    # U = eg_mpc(model,A,B,deepcopy(X),deepcopy(U),ϵ_f)
    U,normdu = eg_mpc_pdip(model,A,B,deepcopy(X),deepcopy(U),ϵ_f)
    @show normdu
    if normdu < 1e-4
        @info "CPAG converged"
        break
    end

end
t_vec = (0:(length(X)-1))*dt*3600


# @infiltrate
# error()
# pme = abs.(eps_hist .- ϵ_f)
# mat"
# figure
# hold on
# %plot($eps_hist)
# plot(0:($T-1),$pme,'linewidth',4)
# set(gca, 'YScale', 'log')
# xlabel('Predictor-corrector iteration number')
# ylabel('Terminal specific energy error')
# %set(gca,'FontSize',15)
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
X1m = Xm[1:8,:]
X2m = Xm[9:16,:]
X3m = Xm[17:24,:]

X1 = vec_from_mat(X1m)
X2 = vec_from_mat(X2m)
X3 = vec_from_mat(X3m)
mat"
figure
hold on
title('Bank Angle')
plot($t_vec, rad2deg($Xm(7,:)),'linewidth',4)
xlabel('Time (s)')
xlim([0,$t_vec(end)])
ylabel('Bank Angle (degrees)')
set(gca,'FontSize',14)
hold off
saveas(gcf,'bank.eps','epsc')
"

alt1 = [norm(X1[i][1:3])*1e3-model.evmodel.planet.R for i = 1:length(X)]
alt2 = [norm(X2[i][1:3])*1e3-model.evmodel.planet.R for i = 1:length(X)]
alt3 = [norm(X3[i][1:3])*1e3-model.evmodel.planet.R for i = 1:length(X)]

mat"
figure
hold on
title('Altitude')
plot($t_vec, $alt1)
plot($t_vec, $alt2)
plot($t_vec, $alt3)
legend('f_1','f_2','f_3')
xlabel('Time (s)')
xlim([0,$t_vec(end)])
ylabel('Altitude (km)')
set(gca,'FontSize',14)
hold off
saveas(gcf,'alt.eps','epsc')
"

ϵ1 = X1m[8,:]
ϵ2 = X2m[8,:]
ϵ3 = X3m[8,:]

mat"
figure
hold on
title('Specific Orbital Energy')
plot($t_vec,$ϵ1)
plot($t_vec,$ϵ2)
plot($t_vec,$ϵ3)
plot($t_vec,$ϵ_f*ones(length($t_vec),1),'r--')
legend('f_1','f_2','f_3','target orbit')
xlim([0,$t_vec(end)])
xlabel('Time (s)')
ylabel('Specific Orbital Energy')
set(gca,'FontSize',14)
hold off
saveas(gcf,'eps.eps','epsc')
"

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
