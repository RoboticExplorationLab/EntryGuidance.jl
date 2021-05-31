using LinearAlgebra, Convex, Mosek, MosekTools, ForwardDiff
using MATLAB
using Attitude

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
    x2 = Ad*x + Bd*u
    xdot = x2-x
    return dot(∇h(x),xdot)
end

function cbf(x,u)
    if norm(p - x[1:2])<5
        u_cvx = Variable(2)
        α = .3
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


# double integrator
# x = [r v]
dt = .1
Ac = Array([zeros(2,2) I(2); zeros(2,4) ])
Bc = Array([zeros(2,2);I(2)])

H = exp(dt*[Ac Bc; zeros(2,6)])
Ad = H[1:4,1:4]
Bd = H[1:4,5:6]

# we are going to drive the system to the origin first
N = 200
X = [zeros(4) for i = 1:N]
X[1] = [10;20;0;0]
K = [1*I(2) 2*I(2)]
for i = 1:(N-1)
    u_vanil = -K*X[i]
    u = cbf(X[i],u_vanil)
    # if norm([7;6] - X[i][1:2])<5
    #     @show hdot(X[i],u)
    # end
    X[i+1] = Ad*X[i] + Bd*u
end

xm = mat_from_vec(X)
p = float(p)
mat"
figure
hold on
plot($xm(1,:),$xm(2,:))
plot($xm(1,:),$xm(2,:),'.','markersize',20)
r = 5;
disp($p)

x = $p(1)
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
