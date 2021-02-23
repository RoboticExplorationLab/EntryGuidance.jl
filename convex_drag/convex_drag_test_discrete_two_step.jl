using Convex, Mosek, MosekTools, LinearAlgebra
const cvx = Convex

# parameters
Cd = 2.2
Area = 4.5
ρ = 1.22e-5
g = zeros(3)
# g = [0;0;-3.71]


# discrete dynamics for double integrator
dt = 0.1
Ad = Array(float(I(6)));Ad[1:3,4:6]=dt*Array(I(3))
Bd = [.5*dt^2*diagm(ones(3));diagm(ones(3))]

# initial conditions
rk = [1000;2000;3000]
vk = [500;500;10]

N = 60
# declare variables for a 2 step simulation
d = cvx.Variable(3,N-1) # drag acceleration history
Γ = cvx.Variable(N-1)   # slack variable
r = cvx.Variable(3,N) # position history
v = cvx.Variable(3,N) # velocity history

# initial condition constraint
cons_list = Constraint[]
push!(cons_list,r[:,1] == rk)
push!(cons_list,v[:,1] == vk)

# two dynamics steps
for i = 1:N-1
    # dynamics constraints
    push!(cons_list, [r[:,i+1];v[:,i+1]] == Ad*[r[:,i];v[:,i]] + Bd*(d[:,i]+g))

    # we're trying to replicate this: |d| ≤ .5*ρ*Cd*A|v|² (non-convex)
    push!(cons_list, norm(d[:,i])<= Γ[i])
    push!(cons_list, .5*ρ*Cd*Area*dot(v[:,i],v[:,i]) <= Γ[i])
    #NOTE: our convex inequalities must be [norm(z)<= ...] or [z'*Q*z <= ...]
end

# this is convex and it solves, but there is nothing stopping Γ from going
# to a huge number and allowing unrealistic drag.

# min velocity squared
# problem = cvx.minimize( sumsquares(vec(v)) + 1e6*sum(Γ), cons_list)

# min energy
problem = cvx.minimize( sumsquares(vec(v)) + 1e6*sum(Γ), cons_list)
# problem = cvx.minimize( .5*sumsquares(v[:,end]) + 1e3*sum(Γ), cons_list)
# problem = cvx.minimize( .5*sumsquares(vec(v)) + abs(sum(g))*sumsquares(r[3,:]) + 1e6*sum(Γ), cons_list)
# problem = cvx.minimize(  .5*sumsquares(v[:,end]) + abs(g[3])*sumsquares(r[3,end]) + 1e6*sum(Γ), cons_list)
# problem = cvx.minimize( p + 1e6*sum(Γ), cons_list)
cvx.solve!(problem, () -> Mosek.Optimizer())                 # NOTE: succeeds


# now we get the true solution
r_tru = zeros(3,N)
v_tru = zeros(3,N)
d_tru = zeros(3,N-1)
r_tru[:,1] = copy(rk)
v_tru[:,1] = copy(vk)
Γ_tru = zeros(N-1)
for i = 1:N-1
    vl = @views v_tru[:,i]
    Γ_tru[i] = .5*ρ*Cd*Area*dot(vl,vl)
    d_tru[:,i] = -.5*ρ*Cd*Area*dot(vl,vl)*normalize(vl)
    next_state = Ad*[r_tru[:,i];v_tru[:,i]] + Bd*(d_tru[:,i] + g)
    r_tru[:,i+1] = next_state[1:3]
    v_tru[:,i+1] = next_state[4:6]
end

r_cvx = r.value
v_cvx = v.value
d_cvx = d.value
using MATLAB
mat"
figure
hold on
title('Position')
plot3($r_cvx(1,:),$r_cvx(2,:),$r_cvx(3,:))
plot3($r_tru(1,:),$r_tru(2,:),$r_tru(3,:),'o')
legend('cvx','true')
view(0,0)
hold off
"

@show r.value[3,end]
Γ_cvx = Γ.value

real_d_cvx = [.5*ρ*Cd*Area*dot(v_cvx[:,i],v_cvx[:,i]) for i = 1:size(v_cvx,2)]

# mat"
# figure
# hold on
# title('Slack (drag norm limit)')
# plot($Γ_cvx)
# plot($Γ_tru,'o')
# legend('cvx','true')
# hold off
# "

vnorm_tru = [norm(v_tru[:,i]) for i = 1:size(v_tru,2)]
vnorm_cvx = [norm(v_cvx[:,i]) for i = 1:size(v_cvx,2)]

# mat"
# figure
# hold on
# title('Velocity norm')
# plot($vnorm_cvx)
# plot($vnorm_tru)
# legend('cvx','true')
# hold off
# "

dnorm_tru = [norm(d_tru[:,i]) for i = 1:size(d_tru,2)]
dnorm_cvx = [norm(d_cvx[:,i]) for i = 1:size(d_cvx,2)]

mat"
figure
hold on
title('Drag Norms from CVX')
plot($real_d_cvx(1:end-1))
plot($Γ_cvx,'o')
plot($dnorm_cvx,'--')
hold off
"

# mat"
# figure
# hold on
# title('drag norm')
# plot($dnorm_cvx)
# plot($dnorm_tru)
# legend('cvx','true')
# hold off
# "

mat"
figure
hold on
title('Velocity')
subplot(3,1,1)
hold on
plot($v_cvx(1,:))
plot($v_tru(1,:))
ylabel('Vx')

subplot(3,1,2)
hold on
plot($v_cvx(2,:))
plot($v_tru(2,:))
ylabel('Vy')
subplot(3,1,3)
hold on
plot($v_cvx(3,:))
plot($v_tru(3,:))
ylabel('Vz')
hold off
"


# get direction of drag and velocity

# normalize(v_cvx[:,1])
# normalize(d_cvx[:,1])


# quiver plot
mat"
figure
hold on
quiver3($r_cvx(1,1:end-1),$r_cvx(2,1:end-1),$r_cvx(3,1:end-1),$d_cvx(1,:),$d_cvx(2,:),$d_cvx(3,:))
quiver3($r_tru(1,1:end-1),$r_tru(2,1:end-1),$r_tru(3,1:end-1),$d_tru(1,:),$d_tru(2,:),$d_tru(3,:))
legend('CVX','Sim')
view(-15,5)
hold off
"
