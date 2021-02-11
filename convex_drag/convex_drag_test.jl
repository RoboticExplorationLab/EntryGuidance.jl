using Convex, Mosek, MosekTools, LinearAlgebra, COSMO

const cvx = Convex

v = randn(3)
Cd = 2.2
A = 4.5
ρ = 1.22e-1

max_drag_norm = .5*Cd*ρ*A*dot(v,v)

real_drag = max_drag_norm*(-normalize(v))

d = cvx.Variable(3)

problem = cvx.minimize(dot(d,v),[norm(d)<=max_drag_norm])

# cvx.solve!(problem, () -> COSMO.Optimizer())
cvx.solve!(problem, () -> Mosek.Optimizer())

@test isapprox(vec(d.value),real_drag,atol = 1e-6)

@show norm(vec(d.value) - real_drag)
