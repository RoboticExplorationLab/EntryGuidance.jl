using LinearAlgebra
using EntryGuidance
const EG = EntryGuidance
using StaticArrays
using Attitude
using MATLAB
using Infiltrator
using ForwardDiff
using Convex
using Mosek, MosekTools
using COSMO
using SuiteSparse
using SparseArrays
using Interpolations


include(joinpath(@__DIR__,"dynamics.jl"))
include(joinpath(@__DIR__,"rollout_stuff.jl"))
include(joinpath(@__DIR__,"mpc.jl"))
include(joinpath(@__DIR__,"post_process.jl"))

function first_test()


alt0= 145.0
γ0 = deg2rad(-15)
V0 = 12.845*3600

# evmodel = CartesianMSLModel()
model = EntryVehicle(CartesianMSLModel(),1.0)
#Initial conditions for MSL
Rm = model.evmodel.planet.R
r0 = [Rm+alt0, 0.0, 0.0] #Atmospheric interface at 125 km altitude
# V0 = 5.845*3600 #Mars-relative velocity at interface of 5.845 km/sec
# γ0 = -15.474*(pi/180.0) #Flight path angle at interface
v0 = V0*[sin(γ0), cos(γ0), 0.0]
epsilon0 = epsilon(model,[r0;v0])
σ0 = deg2rad(30)
x0 = [r0;v0;σ0;epsilon0]
# x0 = [3443.300786841311, 270.4345771068569, 0.0, -6051.64651501579, 20222.23824790719, 0.0]

# first rollout
dt = 2/3600
N = 100
X = NaN*[@SArray zeros(8) for i = 1:N]
U = [@SArray zeros(1) for i = 1:N-1]

X[1] = deepcopy(x0)
end_idx = NaN
for i = 1:(N-1)
    U[i] = [0]
    X[i+1] = rk4(model,X[i],U[i],dt)
end


# let's try some CPAG stuff
# rollout(model,x0,U_in,dt)
# -0.4 is the goal epsilon
ϵ_f = -0.4

T = 10
eps_hist = zeros(T)
for i = 1:T

    X = rollout(model,x0,U,dt)
    eps_hist[i] = X[end][8]

    # linearize
    A,B = getAB(model,X,U,dt)

    # solve cvx prob (correct)
    U = eg_mpc(model,A,B,deepcopy(X),deepcopy(U),ϵ_f)

end



mat"
figure
hold on
plot($eps_hist)
hold off
"

eps2 = zeros(length(X))
for i = 1:length(X)
    eps2[i] = epsilon(model,X[i])
end

mat"
figure
hold on
plot($eps2)
hold off
"

Xm = mat_from_vec(X)
mat"
figure
hold on
title('Bank Angle')
plot(rad2deg($Xm(7,:)))
hold off
"

alt = [norm(X[i][1:3])-model.evmodel.planet.R for i = 1:length(X)]

mat"
figure
hold on
title('Altitude')
plot($alt)
hold off
"

# Um = mat_from_vec(U)
#
# mat"
# figure
# hold on
# plot($Um')
# hold off
# "



# --------- general plotting -----------
# alt = [norm(X[i][1:3])-model.evmodel.planet.R for i = 1:length(X)]
# #
# epsilon2 = zeros(length(X))
# for i = 1:length(X)
#     μ =  model.evmodel.planet.gravity.μ
#     r = X[i][1:3]
#     v = X[i][4:6]
#     epsilon2[i] = (dot(v,v)/2 - μ/norm(r))/1e8
# end
# # # @infiltrate
# #
# mat"
# figure
# hold on
# title('Altitude')
# plot($alt)
# hold off
# "
#
# mat"
# figure
# hold on
# title('post process epsilon')
# plot($epsilon2)
# hold off
# "
#
# Xm = mat_from_vec(X)
# eps = Xm[8,:]
# mat"
# figure
# hold on
# title('dynamics epsilon')
# plot($eps)
# hold off
# "

    return nothing
end
first_test()



# mat"
# figure
# hold on
# plot($xm(1:3,:)')
# hold off
# "
#
# mat"
# figure
# hold on
# plot($tt,$xm(3,:),'o')
# plot($tt,$z2)
# hold off
# "


# aa = [randn(3,10) for i = 1:3]
#
# mat"
# figure
# hold on
# for i = 1:length($aa)
#     p = $aa{1}
#     disp(p)
# end
# hold off
# "
#
