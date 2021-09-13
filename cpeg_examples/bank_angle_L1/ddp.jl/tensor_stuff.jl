using LinearAlgebra, TensorCore, ForwardDiff
const FD = ForwardDiff
# Gxx = Q + A'*P*A + p' ⊡ fxx
# Guu = R + B'*P*B + p' ⊡ fuu
# Gux = B'*P*A     + p' ⊡ fux
model = EntryVehicle(CartesianMSLModel(),1e4)
#Initial conditions for MSL
Rm = model.evmodel.planet.R
r0 = [Rm+125.0, 0.0, 0.0] #Atmospheric interface at 125 km altitude
V0 = 5.845*3600 #Mars-relative velocity at interface of 5.845 km/sec
γ0 = -15.474*(pi/180.0) #Flight path angle at interface
v0 = V0*[sin(γ0), cos(γ0), 0.0]
x0 = [r0;v0]

#Final conditions for MSL 631.979 km down range and 7.869 km cross-range from entry point
Rf = Rm+10.0 #Parachute opens at 10 km altitude
rf = Rf*[cos(7.869/Rf)*cos(631.979/Rf); cos(7.869/Rf)*sin(631.979/Rf); sin(7.869/Rf)]
vf = zeros(3)
xf = [rf;v0/10]


x = copy(xf)
u = 100*randn(2)
dt = 2/3600
A = ForwardDiff.jacobian(dx->rk4(model,dx,u,dt),x)
B = ForwardDiff.jacobian(du->rk4(model,x,du,dt),u)


dAdx = ForwardDiff.jacobian(_x -> Avec(model,_x,u,dt),x)
dBdx = ForwardDiff.jacobian(_x -> Bvec(model,_x,u,dt),x)
@btime dBdu = ForwardDiff.jacobian(_u -> Bvec(model,x,_u,dt),u)

fxx = reshape(dAdx,6,6,6)
fuu = reshape(dBdu,6,2,2)
fux = reshape(dBdx,6,2,6)

P = randn(6,6)
Q = randn(6,6)
q = randn(6)
R = randn(2,2)
p = randn(6)
Rxx = p' ⊡ fxx
Ruu = p' ⊡ fuu
Rux = p' ⊡ fux


Rxx2 = FD.hessian( _x -> dot(p,rk4(model,_x,u,dt)),x)
Ruu2 = FD.hessian( _u -> dot(p,rk4(model,x,_u,dt)),u)
@btime Rux2 = FD.jacobian(_x -> p'*ForwardDiff.jacobian(du->rk4(model,_x,du,dt),u),x)
