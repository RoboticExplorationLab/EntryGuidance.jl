using LinearAlgebra
using EntryGuidance
const EG = EntryGuidance
using StaticArrays
using Attitude
using MATLAB
using Infiltrator
using ForwardDiff
# using Convex
# using Mosek, MosekTools
# using COSMO
using SuiteSparse
using SparseArrays
# using Interpolations


include(joinpath(@__DIR__,"dynamics.jl"))
# include(joinpath(@__DIR__,"rollout_stuff.jl"))
# include(joinpath(@__DIR__,"mpc.jl"))
include(joinpath(@__DIR__,"post_process.jl"))
include(joinpath(@__DIR__,"srekf.jl"))

function first_test()


alt0= 145.0
γ0 = deg2rad(-15)
V0 = 12.845*3600

# evmodel = CartesianMSLModel()
dscale = 1.0e3
tscale = 1.0
uscale = 1.0
model = EntryVehicle(CartesianMSLModel(),dscale,tscale,uscale)
#Initial conditions for MSL
Rm = model.evmodel.planet.R
r0 = [Rm+alt0, 0.0, 0.0] #Atmospheric interface at 125 km altitude


# H =7.295#, 5.25e7
# ρ0 = 5.25
# V0 = 5.845*3600 #Mars-relative velocity at interface of 5.845 km/sec
# γ0 = -15.474*(pi/180.0) #Flight path angle at interface
v0 = V0*[sin(γ0), cos(γ0), 0.0]
σ0 = deg2rad(30)
# H =11.1
# ρ0 = 2.0 # this gets multiplied by 1e7 at one point
kρ = 2.0
x0 = [r0/dscale;v0/(dscale/tscale); σ0;kρ]

# first rollout
dt = 1/3600/tscale
N = 150
X = NaN*[zeros(8) for i = 1:N]
U = [zeros(1) for i = 1:N-1]
Y = [zeros(7) for i = 1:N]

X[1] = deepcopy(x0)
μ = deepcopy(X)
μ[1][8] = μ[1][8] + 0.1*randn()
F = [zeros(8,8) for i = 1:N]
Σ = (0.001*Matrix(float(I(8))))
Σ[8,8] = (1)^2
F[1] = chol(Matrix(Σ))
#t
Q = 1e-20*Matrix(I(8))
# Q[8,8] = 1e-10
# Q = diagm( [(.00005)^2*ones(3)/dscale; .00005^2*ones(3)/(dscale/tscale); (1e-5)^2;(1e-5)^2])
# R = 1e-10*Matrix(I(6))
R = diagm( [(.1)^2*ones(3)/dscale; (0.0002)^2*ones(3)/(dscale/tscale)])

# @show diag(Q)
@show diag(R)
kf_sys = (dt = dt, ΓR = chol(R), ΓQ = chol(Q))

end_idx = NaN
for i = 1:(N-1)
    U[i] = [sin(i/10)/30]
    X[i+1] = rk4(model,X[i],U[i],dt)
    Y[i+1] = measurement(model,X[i+1]) + kf_sys.ΓR*randn(6)

    μ[i+1], F[i+1] = sqrkalman_filter(model, μ[i],F[i],U[i],Y[i+1],kf_sys)
end

Xm = mat_from_vec(X)

mat"
figure
hold on
plot($Xm(7,:))
hold off
"
alt, dr, cr = postprocess(model::EntryVehicle,X,x0)

alt_k, dr_k, cr_k = postprocess(model::EntryVehicle,μ,x0)

perr = ([dscale*norm(X[i][1:3] - μ[i][1:3]) for i = 1:N])
verr = ([(dscale/tscale)*norm(X[i][4:6] - μ[i][4:6]) for i = 1:N])

yperr = ([dscale*norm(X[i][1:3] - Y[i][1:3]) for i = 2:N])
yverr = ([(dscale/tscale)*norm(X[i][4:6] - Y[i][4:6]) for i = 2:N])
# @infiltrate
# error()
# mat"
# figure
# hold on
# plot($alt)
# plot($alt_k)
# hold off
# "
#
# mat"
# figure
# hold on
# plot($dr,$cr)
# plot($dr_k,$cr_k)
# xlabel('downrange')
# ylabel('crossrange')
# hold off
# "
mat"
figure
hold on
title('Position Error')
plot($perr)
plot($yperr)
legend('EKF','Y')
ylabel('km')
hold off
"
mat"
figure
hold on
title('Velocity Error')
plot($verr)
plot($yverr)
legend('EKF','Y')
ylabel('km/s')
hold off
"
μm = mat_from_vec(μ)
mat"
figure
hold on
plot($Xm(8,:)','b')
plot($μm(8,:)','r')
hold off
"

return nothing
end

let
    first_test()
end


#
# alts = 1:150
# ρ1 = zeros(length(alts))
# ρ2 = zeros(length(alts))
# for i = 1:length(alts)
#     ρ1[i] = atmospheric_density([Rm+alts[i], 0.0, 0.0], model.evmodel)
#     ρ2[i] = my_atmo(model,[Rm+alts[i], 0.0, 0.0],H,ρ0)
# end
#
# mat"
# figure
# hold on
# plot($ρ1)
# plot($ρ2)
# legend('BI','Straight')
# set(gca, 'YScale', 'log')
# hold off
# "
