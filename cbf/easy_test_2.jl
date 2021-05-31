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
    # d >= 3 so d - 3 >= 0
    # p = [7;6]
    d = norm(x[1:2] - p)
    radius = 3
    return d - 3
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
    return [zeros(2,2) I(2); zeros(2,4) ]*x + [zeros(2,2);I(2)]*u
end

# we are going to drive the system to the origin first
dt = .1
N = 200
X = [zeros(4) for i = 1:N]
X[1] = [10;20;0;0]
K = [1*I(2) 2*I(2)]
for i = 1:(N-1)
    u_vanil = -K*X[i]
    u = cbf(X[i],u_vanil)
    X[i+1] = rk4(dynamics,X[i],u,dt)
end

xm = mat_from_vec(X)
p = float(p)
mat"
figure
hold on
plot($xm(1,:),$xm(2,:))
plot($xm(1,:),$xm(2,:),'.','markersize',20)
r = 5;

x = $p(1);
y = $p(2);
th = 0:pi/50:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
plot(xunit, yunit)

r = 3;
th = 0:pi/50:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
plot(xunit, yunit)
axis equal
hold off
"
