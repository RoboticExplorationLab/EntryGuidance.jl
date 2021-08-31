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
x0 = [r0;v0]
# x0 = [3443.300786841311, 270.4345771068569, 0.0, -6051.64651501579, 20222.23824790719, 0.0]

# first rollout
dt = 2/3600
N = 100
X = NaN*[@SArray zeros(6) for i = 1:N]
U = [@SArray zeros(2) for i = 1:N-1]

X[1] = deepcopy(x0)
end_idx = NaN
for i = 1:(N-1)
    U[i] = [0;0.5]
    X[i+1] = rk4(model,X[i],U[i],dt)
    if altitude(model,X[i+1])<10
        X = X[1:i]
        @info "under altitude on first rollout"
        break
    end
end

alt = [norm(X[i][1:3])-model.evmodel.planet.R for i = 1:length(X)]

epsilon = zeros(length(X))
for i = 1:length(X)
    μ =  model.evmodel.planet.gravity.μ
    r = X[i][1:3]
    v = X[i][4:6]
    epsilon[i] = dot(v,v)/2 - μ/norm(r)
end
# @infiltrate

mat"
figure
hold on
plot($alt)
hold off
"

mat"
figure
hold on
plot($epsilon)
hold off
"

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
