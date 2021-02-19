using Convex, Mosek, MosekTools, LinearAlgebra, COSMO
const cvx = Convex

# parameters
Cd = 2.2
Area = 4.5
ρ = 1.22e-4

# discrete dynamics for double integrator
Ad = Array(float(I(6)));Ad[1:3,4:6]=dt*Array(I(3))
Bd = [.5*dt^2*diagm(ones(3));diagm(ones(3))]

# initial conditions
rk = 100*randn(3)
vk = 100*randn(3)

# declare variables for a 2 step simulation
d = cvx.Variable(3,2) # drag acceleration history
Γ = cvx.Variable(2)   # slack variable
r = cvx.Variable(3,3) # position history
v = cvx.Variable(3,3) # velocity history

# initial condition constraint
cons_list = Constraint[]
push!(cons_list,r[:,1] == rk)
push!(cons_list,v[:,1] == vk)

# two dynamics steps
for i = 1:2
    # dynamics constraints
    push!(cons_list, [r[:,i+1];v[:,i+1]] == Ad*[r[:,i];v[:,i]] + Bd*d[:,i])

    # we're trying to replicate this: |d| ≤ .5*ρ*Cd*A|v|² (non-convex)
    push!(cons_list, norm(d[:,i])<= Γ[i])
    push!(cons_list, .5*ρ*Cd*Area*dot(v[:,i],v[:,i]) <= Γ[i])
    #NOTE: our convex inequalities must be [norm(z)<= ...] or [z'*Q*z <= ...]
end

# this is convex and it solves, but there is nothing stopping Γ from going
# to a huge number and allowing unrealistic drag.
problem = cvx.minimize( sumsquares(vec(v)), cons_list)
cvx.solve!(problem, () -> COSMO.Optimizer(max_iter = 100000))
