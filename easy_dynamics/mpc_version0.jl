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

x0 = [100,200,1000,13,15,-20.0]
Nr = 100
Uc = [zeros(2) for i = 1:2*Nr]
dt = 0.5
X, t_vec, t_impact = rollout(x0,Uc,dt)
# U = U[1:(length(X)-1)]

# xm = mat_from_vec(tX)
T = 70
Xsim = [zeros(6) for i = 1:T]
Xsim[1] = x0
Usim = [zeros(2) for i = 1:T-1]

    # MPC loop
    for i = 1:T-1

        # rollout current plan and find out when we dying
        Xr, Ur, t_vec, t_impact = rollout(Xsim[i],Uc,dt)

        @show length(Xr)
        @show length(Ur)
        @show Xr[end][3]

        # jacobians
        A,B = getAB(Xr,Ur,dt)

        # MPC solve
        Xc, Uc = eg_mpc(A,B,Xr,Ur)

        # @infiltrate
        # testing stuff
        cvxX = mat_from_vec(Xc)
        cvxU = mat_from_vec(Uc)
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
        # @show cvxX[3,end]
        # @infiltrate
        # if i == 10
        #     error()
        # end

        # actual dynamics
        Usim[i] = copy(Uc[1])

        Xsim[i+1] = rk4(Xsim[i],Usim[i],dt)

    end
    xm = mat_from_vec(Xsim)
    mat"
    figure
    hold on
    title('Positions')
    plot($xm(1:3,:)')
    hold off
    "
    mat"
    figure
    hold on
    title('Positions')
    plot($xm(1,:),$xm(2,:))
    plot($xm(1,1),$xm(2,1),'r*')
    hold off
    "
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
