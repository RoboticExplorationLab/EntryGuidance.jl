using LinearAlgebra, Convex, Mosek, MosekTools, ForwardDiff
using MATLAB
using Attitude
const FD = ForwardDiff
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
# now let's add an obstacle (circle) at (4,6) of radius 3
p = [6;6]
function h_safety(x)
    θ = 60
    avoid_N = normalize([cosd(θ);sind(θ);0])

    camera_B = [1;0;0]

    N_Q_B = dcm_from_p( x[1:3])
    camera_N = (N_Q_B)*camera_B

    return acos(dot(avoid_N,camera_N)) - deg2rad(30)
end

function ∇h(x)
    return ForwardDiff.gradient(h_safety,x)
end
function hdot(x,u)
    nu = size(u,1)
    A,B = jacobians(dynamics,x,zeros(nu),dt)
    x2 = rk4(dynamics,x,zeros(nu),dt) + B*u
    xdot = x2-x
    return dot(∇h(x),xdot)
end

function cbf(x,u)
    if h_safety(x)<2.0
        u_cvx = Variable(length(u))
        α = .2
        problem = minimize(sumsquares(u-u_cvx), [hdot(x,u_cvx) >= -α*h_safety(x)])
        solve!(problem, () -> Mosek.Optimizer())
        return evaluate(u_cvx)
    else
        return u
    end
end
function dynamics(x,u)
    p = x[1:3]
    ω = x[4:6]
    # qdot = .5*(q ⊙ [0;ω])
    pdot = pdot_from_w(p,ω)
    J = Diagonal([1;2;3])
    ω̇ = J\(u - (ω × (J*ω)))
    return [pdot;ω̇]
end

# we are going to drive the system to the origin first
dt = .1
N = 200
X = [zeros(6) for i = 1:N]
Random.seed!(4)
X[1] = [p_from_phi(deg2rad(140)*normalize(randn(3))); zeros(3)]
h = zeros(N-1)
kp = 1
kd = 5
for i = 1:(N-1)
    # u_vanil = -K*X[i]
    ϕ = phi_from_p(X[i][1:4])
    ω = X[i][4:6]

    u = -(kp*ϕ + kd*ω)
    u = cbf(X[i],u)
    h[i]=h_safety(X[i])

    X[i+1] = rk4(dynamics,X[i],u,dt)
end

xm = mat_from_vec(X)

# phim = mat_from_vec([phi_from_p(X[i][1:3]) for i = 1:N])
phim = (xm[1:3,:])
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
mat"
figure
hold on
title('Safety')
plot($h)
hold off
"
