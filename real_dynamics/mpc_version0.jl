using LinearAlgebra
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


function first_test()

# x0 = [100,200,1000,13,15,-20.0]
x0 = [-1200,-1200,1000,90,150,-20.0]
Nr = 100
Uc = [.1*zeros(2) for i = 1:2*Nr]
dt = 0.5
X, t_vec, t_impact = rollout(x0,Uc,dt)
# U = U[1:(length(X)-1)]

# xm = mat_from_vec(tX)
traj_hist = Array{Array{Float64,1},1}[]
T = 40
Xsim = [zeros(6) for i = 1:T]
Xsim[1] = x0
Usim = [zeros(2) for i = 1:T-1]

    # MPC loop
    for i = 1:T-1

        # rollout current plan and find out when we dying
        Xr, Ur, t_vec, t_impact = rollout(deepcopy(Xsim[i]),Uc[2:end],dt)
        push!(traj_hist,Xr)
        @show length(Xr)
        @show length(Ur)
        @show Xr[end][3]
        @assert Xsim[i] == Xr[1]

        # jacobians
        # Ur = Ur .*0
        A,B = getAB(Xr,Ur,dt)

        # MPC solve
        @show length(Uc)
        @show length(Ur)
        Xc, Uc = eg_mpc(A,B,Xr,Ur)

        # @infiltrate
        # testing stuff
        cvxX = mat_from_vec(Xc)
        cvxU = mat_from_vec(Uc)
        # if i == T-1
        #     @info "end game"
        #     Xr, Ur, t_vec, t_impact = rollout(Xsim[i],Uc,dt)
        #     push!(traj_hist,Xr)
        #     # mat"
        #     # figure
        #     # hold on
        #     # title('Positions')
        #     # plot($cvxX(1:3,:)')
        #     # hold off
        #     # "
        #     # mat"
        #     # figure
        #     # hold on
        #     # title('Velocities')
        #     # plot($cvxX(4:6,:)')
        #     # hold off
        #     # "
        #     # mat"
        #     # figure
        #     # hold on
        #     # title('Controls')
        #     # plot($cvxU')
        #     # hold off
        #     # "
        # end

        # actual dynamics
        Usim[i] = copy(Uc[1])

        Xsim[i+1] = rk4(Xsim[i],Usim[i],dt)

    end
    xm = mat_from_vec(Xsim)
    # mat"
    # figure
    # hold on
    # title('Positions')
    # plot($xm(1:3,:)')
    # hold off
    # "
    # mat"
    # figure
    # hold on
    # title('Positions')
    # plot($xm(1,:),$xm(2,:))
    # plot($xm(1,1),$xm(2,1),'r*')
    # hold off
    # "
    um = mat_from_vec(Usim)
    mat"
    figure
    hold on
    title('Controls')
    plot($um')
    hold off
    "
    # mat"
    # figure
    # hold on
    # title('Positions')
    # plot($cvxX(1:3,:)')
    # hold off
    # "
    # mat"
    # figure
    # hold on
    # title('Velocities')
    # plot($cvxX(4:6,:)')
    # hold off
    # "
    # mat"
    # figure
    # hold on
    # title('Controls')
    # plot($cvxU')
    # hold off
    # "
    # mat"
    # figure
    # hold on
    # title('Controls')
    # plot($cvxU(1,:),$cvxU(2,:))
    # hold off
    # "
    @show length(traj_hist)
    traj = [mat_from_vec(traj_hist[i]) for i = 1:length(traj_hist)]
    # @infiltrate
    # error()
    mat"
    figure
    hold on
    rgb1 = [29 38 113]/255;
    rgb2 = [195 55 100]/255;
    drgb = rgb2-rgb1
    for i = 1:length($traj)
        p = $traj{i};
        plot3(p(1,:),p(2,:),p(3,:),'Color',rgb1 + drgb*(i-1)/length($traj),'linewidth',3)
        plot3(p(1,1),p(2,1),p(3,1),'r.','markersize',30)
    end
    plot3([0],[0],[0],'g.','markersize',30)
    zlim([0,1100])
    %axis equal
    % grid on
    hold off
    "

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
