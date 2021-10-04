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
using COSMO
using SuiteSparse
using SparseArrays
using Interpolations

include(joinpath(@__DIR__,"dynamics.jl"))
include(joinpath(@__DIR__,"rollout_stuff.jl"))
include(joinpath(@__DIR__,"mpc.jl"))
include(joinpath(@__DIR__,"post_process.jl"))

function first_test()

# evmodel = CartesianMSLModel()
model = EntryVehicle(CartesianMSLModel())
#Initial conditions for MSL
Rm = model.evmodel.planet.R
r0 = [Rm+125.0, 0.0, 0.0] #Atmospheric interface at 125 km altitude
V0 = 5.845*3600 #Mars-relative velocity at interface of 5.845 km/sec
γ0 = -15.474*(pi/180.0) #Flight path angle at interface
v0 = V0*[sin(γ0), cos(γ0), 0.0]
σ = deg2rad(0)
# x0 = [r0;v0;σ]

x0 = [3443.300786841311, 270.4345771068569, 0.0, -6051.64651501579, 20222.23824790719, 0.0,σ]

#Final conditions for MSL 631.979 km down range and 7.869 km cross-range from entry point
Rf = Rm+10.0 #Parachute opens at 10 km altitude
rf = Rf*[cos(7.869/Rf)*cos(631.979/Rf); cos(7.869/Rf)*sin(631.979/Rf); sin(7.869/Rf)]
vf = zeros(3)
xf = [rf;vf;0]

# first rollout
dt = 2/3600
N = 180
X = NaN*[@SArray zeros(7) for i = 1:N]
U = [0 for i = 1:N-1]

X[1] = deepcopy(x0)
end_idx = NaN
for i = 1:(N-1)
    X[i+1] = rk4(model,X[i],U[i],dt)
    if altitude(model,X[i+1])<10
        @info "under altitude on first rollout"
        end_idx = i + 1
        break
    end
end
X = X[1:end_idx]
U = U[1:end_idx]
Uc = deepcopy(U)
# @infiltrate
# error()
T = 50
Xsim = [zeros(6) for i = 1:T]
Xsim[1] = x0
Usim = [0.0 for i = 1:T-1]
althist = [zeros(2) for i = 1:T-1]
drhist = [zeros(2) for i = 1:T-1]
crhist = [zeros(2) for i = 1:T-1]
# HRhist = zeros(T-1)
# MPC loop
for i = 1:T-1

    # rollout current plan and find out when we dying
    Xr, Ur, t_vec, t_impact = rollout(model,deepcopy(Xsim[i]),Uc[2:end],dt)
    @assert Xsim[i] == Xr[1]

    althist[i], drhist[i], crhist[i] = postprocess(model::EntryVehicle,Xr,x0)
    # HRhist[i] = get_heating(model::EntryVehicle, Xsim[i])
    # push!(althist,alt)
    # push!(drhist,dr)
    # push!(crhist,cr)
    # jacobians
    A,B = getAB(model,Xr,Ur,dt)

    # MPC solve
    Xc, Uc = eg_mpc(model::EntryVehicle,A,B,Xr,Ur,xf)

    cvxX = mat_from_vec(Xc)
    cvxU = deepcopy(Uc)

    # actual dynamics
    Usim[i] = copy(Uc[1])

    # @infiltrate
    # error()
    Xsim[i+1] = rk4(model::EntryVehicle,Xsim[i],Usim[i],dt)

end
    # xm = mat_from_vec(Xsim)
    # # mat"
    # # figure
    # # hold on
    # # title('Positions')
    # # plot($xm(1:3,:)')
    # # hold off
    # # "
    # # mat"
    # # figure
    # # hold on
    # # title('Positions')
    # # plot($xm(1,:),$xm(2,:))
    # # plot($xm(1,1),$xm(2,1),'r*')
    # # hold off
    # # "
    # um = mat_from_vec(Usim)
    # um = copy(Usim)
    # mat"
    # figure
    # hold on
    # title('Controls')
    # plot($um')
    # hold off
    # "
    xm = copy(mat_from_vec(Xsim))
    # @infiltrate
    bank = rad2deg.(xm[7,:])
    mat"
    figure
    hold on
    %title('Bank Angle')
    plot(0:2:(2*(length($bank)-1)), $bank,'linewidth',4)
    ylabel('bank angle (degrees)')
    xlabel('time (seconds)')
    hold off
    addpath('/Users/kevintracy/devel/WiggleSat/matlab2tikz-master/src')
    %matlab2tikz('bank_angle.tex')
    %saveas(gcf,'bank.png')
    "
    # # mat"
    # # figure
    # # hold on
    # # title('Positions')
    # # plot($cvxX(1:3,:)')
    # # hold off
    # # "
    # # mat"
    # # figure
    # # hold on
    # # title('Velocities')
    # # plot($cvxX(4:6,:)')
    # # hold off
    # # "
    # # mat"
    # # figure
    # # hold on
    # # title('Controls')
    # # plot($cvxU')
    # # hold off
    # # "
    # # mat"
    # # figure
    # # hold on
    # # title('Controls')
    # # plot($cvxU(1,:),$cvxU(2,:))
    # # hold off
    # # "
    # @show length(traj_hist)
    # traj = [mat_from_vec(traj_hist[i]) for i = 1:length(traj_hist)]
    # # @infiltrate
    # # error()
    # althist = [mat_from_vec]
    xf_dr, xf_cr = rangedistances(model,xf,x0)
    mat"
    figure
    hold on
    rgb1 = [29 38 113]/255;
    rgb2 = 1.3*[195 55 100]/255;
    drgb = rgb2-rgb1;
    for i = 1:length($drhist)
        px = $drhist{i};
        py = $crhist{i};
        if i < 8
            plot(px,py,'Color',rgb1 + drgb*(i-1)/7,'linewidth',3)
        end
        %plot(px(1),py(1),'r.','markersize',20)
    end
    plot($xf_dr,$xf_cr,'g.','markersize',20)
    xlabel('downrange (km)')
    ylabel('crossrange (km)')
    hold off
    %saveas(gcf,'range.png')
    addpath('/Users/kevintracy/devel/WiggleSat/matlab2tikz-master/src')
    %matlab2tikz('bank_track.tex')
    %close all
    "
    mat"
    figure
    hold on
    rgb1 = [29 38 113]/255;
    rgb2 = 1.3*[195 55 100]/255;
    drgb = rgb2-rgb1;
    for i = 1:length($althist)
        px = $drhist{i};
        alt = $althist{i};
        if i < 6
            colo = drgb*(i-1)/5
            plot(px,alt,'Color',rgb1 + colo,'linewidth',3)
        end
        %plot(px(1),alt(1),'r.','markersize',20)
    end
    plot([0,800],ones( 2,1)*10,'r' )
    plot($xf_dr,10,'g.','markersize',20)
    xlim([100 450])
    xlabel('downrange (km)')
    ylabel('altitude (km)')
    hold off
    %saveas(gcf,'alt.png')
    addpath('/Users/kevintracy/devel/WiggleSat/matlab2tikz-master/src')
    %matlab2tikz('bank_alt.tex')
    "

    # mat"
    # figure
    # plot($HRhist)
    # hold on
    # "

    # AoA, bank = processU(model::EntryVehicle,Xsim,Usim)
    # mat"
    # figure
    # hold on
    # title('Angle of Attack')
    # plot(rad2deg($AoA))
    # hold off
    # "
    # mat"
    # figure
    # hold on
    # title('Bank Angle')
    # plot(rad2deg($bank))
    # hold off
    # "
    return Xsim
end

Xsim2 = first_test()



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
