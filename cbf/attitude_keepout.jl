using LinearAlgebra, Convex, Mosek, MosekTools, ForwardDiff
using MATLAB
using Attitude
const FD = ForwardDiff
using Random
function dynamics(x,u)
    q = x[1:4]
    ω = x[5:7]
    qdot = .5*(q ⊙ [0;ω])
    J = Diagonal([1;2;3])
    ω̇ = J\(u - (ω × (J*ω)))
    return [qdot;ω̇]
end

function rk4(f, x_n, u, h)
    k1 = h*f(x_n,u)
    k2 = h*f(x_n+k1/2, u)
    k3 = h*f(x_n+k2/2, u)
    k4 = h*f(x_n+k3, u)
    return (x_n + (1/6)*(k1+2*k2+2*k3 + k4))
end

function jacobians(f::Function, x, u, dt)
    A = FD.jacobian( _x -> rk4(dynamics, _x,  u, dt), x)
    B = FD.jacobian( _u -> rk4(dynamics,  x, _u, dt), u)
    return A, B
end

function keepout(x)
    θ = 60
    avoid_N = normalize([cosd(θ);sind(θ);0])

    camera_B = [1;0;0]

    N_Q_B = dcm_from_q( x[1:4])
    camera_N = (N_Q_B)*camera_B

    return acos(dot(avoid_N,camera_N)) - deg2rad(40)
end

function runme()


    N = 300
    dt = .1
    X = [zeros(7) for i = 1:N]
    # Random.seed!(6)
    @info "start"
    X[1] = [q_from_phi(deg2rad(140)*[0;0;1]); zeros(3)]

    kp = 1
    kd = 5

    for i = 1:(N-1)

        ϕ = phi_from_q(X[i][1:4])
        ω = X[i][5:7]

        u = -(kp*ϕ + kd*ω)
        X[i+1] = rk4(dynamics,X[i],u,dt)

        d = keepout(X[i])

        if d < 0
            @show i
            @show d
        end

    end

    xm = mat_from_vec(X)

    phim = mat_from_vec([phi_from_q(X[i][1:4]) for i = 1:N])
    t = 0:dt:(N-1)*dt
    # mat"
    # figure
    # title('Axis Angle Attitude')
    # hold on
    # plot($t,rad2deg($phim'))
    # hold off
    # "
    # mat"
    # figure
    # title('Angular Velocity')
    # hold on
    # plot($t,rad2deg($xm(5:7,:)'))
    # hold off
    # "
end

runme()
