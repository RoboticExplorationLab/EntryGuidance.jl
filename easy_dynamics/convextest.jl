using Convex, JuMP, COSMO, Mosek, MosekTools
using JLD2
using Test

function clamp_norm(x,x_max)
    return norm(x)<x_max ? x : normalize(x)*x_max
end
function get_max_L(x)
    v = x[4:6]
    Area = 1.0  # m²
    m = 1000 # kg
    ρ = 1.2*10 # kg/m³
    # Cl = 1.42*α (in radians)
    Cl = 1.42*deg2rad(30)
    # this is maximum allowable lift
    return 0.5*Cl*ρ*Area*dot(v,v)/m
end

# @load "mpc.jld2" mpc

# A,B,X,U = mpc.A, mpc.B, mpc.X, mpc.U
x0 = [-1200,-1200,1000,130,150,-20.0]
dt = 0.5
# N = length(X)
# x0 = copy(X[1])

# first we confirm the rollout
newX = [zeros(6) for i = 1:length(X)]
newU = [zeros(2) for i = 1:length(X)-1]
newX[1] = deepcopy(X[1])
for i = 1:(length(X)-1)
    newU[i] = clamp_norm(U[i],get_max_L(newX[i]))
    newX[i+1] = rk4(newX[i],newU[i],dt)
end

@test norm(X-newX) ≈ 0
@test norm(U-newU) ≈ 0


# max lift vector for each time
L_max = zeros(N)
for i = 1:N
    L_max[i] = get_max_L(X[i])
end

# starting δx
δx0 = x0 - X[1]
@assert norm(x0 - X[1]) ≈ 0
# variables
δx = Variable(6,N)
δu = Variable(2,N-1)

# dynamics constraints
cons = Constraint[ δx[:,1]==δx0 ]

for i = 1:N-1
    push!(cons, δx[:,i+1] == sparse(A[i])*δx[:,i] + B[i]*δu[:,i])
end

# lift vector constraints
for i = 1:N-1
    push!(cons, norm(U[i] + δu[:,i]) <= L_max[i])
    # push!(cons, sumsquares(U[i] + δu[:,i]) <= L_max[i]^2)
end

# problem = minimize(norm(X[N][1:2] + δx[1:2,N]), cons)
α = 1e1
β = 1e2
# problem = minimize(α*norm(X[N][1:2] + δx[1:2,N]) + norm(vec(mat_from_vec(U) + δu)) + norm(vec(δu)), cons)
problem = minimize(α*norm(X[N][1:2] + δx[1:2,N]) + norm(vec(mat_from_vec(U) + δu)), cons)
Convex.solve!(problem, () -> COSMO.Optimizer(eps_abs = 1e-12))

cX = vec_from_mat(mat_from_vec(X) + evaluate(δx))
cU = vec_from_mat(mat_from_vec(U) + evaluate(δu))

dx = δx.value
du = δu.value

mat"
figure
hold on
title('Delta Position')
plot($dx(1:3,:)')
hold off
"
mat"
figure
hold on
title('Delta Velocity')
plot($dx(4:6,:)')
hold off
"
mat"
figure
hold on
title('Delta Control')
plot($du')
hold off
"

cXm = mat_from_vec(cX)
cXu = mat_from_vec(cU)
mat"
figure
hold on
title('Solved Full Position')
plot($cXm(1:3,:)')
hold off
"
mat"
figure
hold on
title('Solved Full Control')
plot($cXu')
hold off
"
