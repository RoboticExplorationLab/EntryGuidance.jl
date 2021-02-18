using Convex, Mosek, MosekTools, LinearAlgebra, COSMO

const cvx = Convex

# v = randn(3)
Cd = 2.2
Area = 4.5
ρ = 1.22e-4

A = [zeros(3,3) Array(I(3));zeros(3,6)]
Bd = [zeros(3,3);Array(I(3))]
dt = 1
Hd = exp([A Bd; zeros(3,9)]*dt)
Ad = Hd[1:6,1:6]
Bd = Hd[1:6,7:9]

rk = randn(3)
vk = randn(3)

# optimization part
max_drag_norm = .5*Cd*ρ*Area*dot(vk,vk)


d = cvx.Variable(3)
Γ = cvx.Variable(1)
vkp1 = cvx.Variable(3)

# dynamics constraint
dyn_con = vkp1 == Ad[4:6,:]*[rk;vk] + Bd[4:6,:]*d

α = .5*ρ*Cd*Area
# drag constraints
con1 = norm(d)          <= Γ
con2 = α*dot(vkp1,vkp1) <= Γ
problem = cvx.minimize(dot(vkp1,vkp1),[con1,
                                       con2,
                                       dyn_con])

# cvx.solve!(problem, () -> COSMO.Optimizer())
cvx.solve!(problem, () -> Mosek.Optimizer())

# @test isapprox(vec(d.value),real_drag,atol = 1e-6)
#
real_drag = max_drag_norm*(-normalize(vk))
@show norm(vec(d.value) - real_drag)/norm(real_drag)
