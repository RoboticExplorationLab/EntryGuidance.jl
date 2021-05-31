using LinearAlgebra, Convex, Mosek, MosekTools, ForwardDiff
using MATLAB
using Attitude
const FD = ForwardDiff
using Random
function dynamics(x,u)
    p = x[1:3]
    ω = x[4:6]
    # qdot = .5*(q ⊙ [0;ω])
    pdot = pdot_from_w(p,ω)
    J = Diagonal([1;2;3])
    ω̇ = J\(u - (ω × (J*ω)))
    return [pdot;ω̇]
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
function h_safety(x)
    θ = 60
    avoid_N = normalize([cosd(θ);sind(θ);0])

    camera_B = [1;0;0]

    N_Q_B = dcm_from_p( x[1:3])
    camera_N = (N_Q_B)*camera_B

    return acos(dot(avoid_N,camera_N)) - deg2rad(40)
end

function ∇h(x)
    return ForwardDiff.gradient(h_safety,x)
end
function hdot(x,u)
    # x2 = Ad*x + Bd*u
    dt = .1
    A,B = jacobians(dynamics,x,zeros(3),dt)
    Δx =  rk4(dynamics,x,zeros(3),dt) + B*u
    return dot(∇h(x),Δx)
end

function cbf(x,u)
    if h_safety(x)<.3
        u_cvx = Variable(3)
        α = .1
        problem = minimize(norm(u-u_cvx), [hdot(x,u_cvx) >= -α*h_safety(x)])
        # problem = minimize(sumsquares(u-u_cvx), [hdot(x,u_cvx) >= 0])
        solve!(problem, () -> COSMO.Optimizer())
        # @infiltrate

        # error()
        return evaluate(u_cvx)
    else
        return u
    end
end
function keepout(x)
    θ = 60
    avoid_N = normalize([cosd(θ);sind(θ);0])

    camera_B = [1;0;0]

    N_Q_B = dcm_from_p( x[1:3])
    camera_N = (N_Q_B)*camera_B

    return acos(dot(avoid_N,camera_N)) - deg2rad(40)
end

function runme()


    N = 300
    dt = .1
    X = [zeros(6) for i = 1:N]
    # Random.seed!(6)
    @info "start"
    X[1] = [p_from_phi(deg2rad(140)*[0;0;1]); zeros(3)]

    kp = 1
    kd = 5

    for i = 1:(N-1)

        ϕ = phi_from_p(X[i][1:4])
        ω = X[i][4:6]

        u = -(kp*ϕ + kd*ω)
        u = cbf(deepcopy(X[i]),u)
        X[i+1] = rk4(dynamics,X[i],u,dt)

        d = keepout(X[i])

        if d < 0
            @show i
            @show d
        end

    end

    xm = mat_from_vec(X)

    phim = mat_from_vec([phi_from_p(X[i][1:3]) for i = 1:N])
    t = 0:dt:(N-1)*dt
    mat"
    figure
    title('Axis Angle Attitude')
    hold on
    plot($t,rad2deg($phim'))
    legend('phi_x','phi_y','phi_z')
    hold off
    "
    mat"
    figure
    title('Angular Velocity')
    hold on
    plot($t,rad2deg($xm(4:6,:)'))
    legend('omega_x','omega_y','omega_z')
    hold off
    "
end

runme()
