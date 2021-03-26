using LinearAlgebra, ForwardDiff
using EntryGuidance
const EG = EntryGuidance
using MATLAB
using Attitude
using StaticArrays
using JuMP, Mosek, MosekTools
using COSMO


struct EntryVehicle{T}
    evmodel::EG.CartesianModel{T}
    uscale::Float64
end

function postprocess(model::EntryVehicle,X,x0)
    N = length(X)
    alt = zeros(N)
    cr = zeros(N)
    dr = zeros(N)
    for i = 1:N
        alt[i] = altitude(model,X[i])
        dr[i],cr[i] = rangedistances(model,X[i],x0)
    end
    return alt, dr, cr
end
function altitude(model,x)
    return norm(x[1:3]) - model.evmodel.planet.R
end
function anglebetween(r1,r2)
    dp = dot(normalize(r1),normalize(r2))
    if dp >1.000000001
        error("over 1 in angle between")
    end
    if dp>1
        dp = 1
    end
    return acos(dp)
end
function rangedistances(model,x,x0)

    # first we get the angular momentum
    r0 = x0[1:3]
    v0 = x0[4:6]
    h = normalize(cross(r0,v0))
    R = model.evmodel.planet.R
    r = x[1:3]

    # downrange stuff
    r_dr = r - dot(r,h)*h
    θ_dr = anglebetween(r0,r_dr)
    dr = θ_dr*R

    # cross range stuff
    dr_r_r = r - r_dr
    θ_cr = anglebetween(r_dr,r)
    cr = dot(dr_r_r,h) > 0 ? θ_cr*R : -θ_cr*R
    return dr, cr

end

# evmodel = CartesianMSLModel()
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
xf = [rf;vf]

function getmaxL(model,x)
    r = x[1:3]
    v = x[4:6]
    ρ = atmospheric_density(r, model.evmodel)
    A = model.evmodel.vehicle.A
    m = model.evmodel.vehicle.m
    Cl = lift_coefficient(deg2rad(30), model.evmodel)
    L = 0.5*Cl*ρ*A*dot(v,v)/m
    return L/model.uscale
end

function evdynamics(model::EntryVehicle, x, u)

    #unpack state
    # r = x[1:3]
    # v = x[4:6]
    # r = SA[x[1],x[2],x[3]]
    # v = SA[x[4],x[5],x[6]]
    r = x[SA[1,2,3]]
    v = x[SA[4,5,6]]

    #atmospheric density
    ρ = atmospheric_density(r, model.evmodel)

    #Calculate drag acceleration
    # Cd = drag_coefficient(α, model.evmodel)
    Cd = drag_coefficient(deg2rad(10), model.evmodel)
    A = model.evmodel.vehicle.A
    m = model.evmodel.vehicle.m
    D = 0.5*Cd*ρ*A*dot(v,v)/m

    #Calculate lift acceleration
    # Cl = lift_coefficient(α, model.evmodel)
    # L = 0.5*Cl*ρ*A*V*V/m

    #get gravity
    g = gravitational_acceleration(r, model.evmodel)

    #Terms involving Ω (planet rotation) are often thrown out in the literature.
    # Ω = model.evmodel.planet.Ω #Set Ω = 0.0 here if you want that behavior
    # Ω̂ = hat([0, 0, Ω])

    #Aerodynamic acceleration
    e1 = cross(r,v)
    e1 = e1/norm(e1)
    e2 = cross(v,e1)
    e2 = e2/norm(e2)
    D_a = -(D/norm(v))*v #+ L*sin(σ)*e1 + L*cos(σ)*e2
    L_a = e1*u[1] + e2*u[2]
                      # this is rotating planet effects
    v̇ = D_a + model.uscale*L_a + g #- 2*Ω̂*v - Ω̂*Ω̂*r

    # return [v; v̇]
    return SA[v[1],v[2],v[3],v̇[1],v̇[2],v̇[3]]
end

function rk4(model,x_n,u,dt)
    k1 = dt*evdynamics(model,x_n,u)
    k2 = dt*evdynamics(model,x_n+k1/2,u)
    k3 = dt*evdynamics(model,x_n+k2/2,u)
    k4 = dt*evdynamics(model,x_n+k3,u)
    return (x_n + (1/6)*(k1 + 2*k2 + 2*k3 + k4))
end
function getAB(model,X,U,dt)
    N = length(X)
    n = length(X[1])
    m = length(U[1])
    A = [@SArray zeros(n,n) for i = 1:(N-1)]
    B = [@SArray zeros(n,m) for i = 1:(N-1)]
    for k = 1:(N-1)
        A[k] = ForwardDiff.jacobian(_x -> rk4(model,_x,U[k],dt),X[k])
        B[k] = ForwardDiff.jacobian(_u -> rk4(model,X[k],_u,dt),U[k])
    end
    return A,B
end

dt = 2/3600
N = 160
X = NaN*[@SArray zeros(6) for i = 1:N]
U = [@SArray zeros(2) for i = 1:N-1]

X[1] = deepcopy(x0)
# global end_idx = NaN
for i = 1:(N-1)
    U[i] = getmaxL(model,X[i])*[0;.5]
    X[i+1] = rk4(model,X[i],U[i],dt)
    if altitude(model,X[i+1])<10
        @info "under altitude"
        # end_idx = i +1
        # X = X[1:(i+1)];U = U[1:(i+1)];
        break
    end
end

# trim
# X = X[1:end_idx];U = U[1:end_idx];
# A,B = getAB(model,X,U,dt)
# quick post process
# alt = [altitude(model,X[i]) for i = 1:length(X)]
xm = mat_from_vec(X)
alt, dr, cr = postprocess(model::EntryVehicle,X,x0)

# mat"
# figure
# hold on
# plot($xm(1:3,:)')
# hold off
# "

mat"
figure
hold on
plot($alt)
plot(1:length($alt),ones( length($alt),1)*10,'r' )
hold off
"


xf_dr, xf_cr = rangedistances(model,xf,x0)
mat"
figure
hold on
plot($dr,$cr)
plot($xf_dr,$xf_cr,'r.','markersize',20)
xlabel('Downrange')
ylabel('Crossrange')
legend('Trajectory','Goal')
hold off
"
#
# um = mat_from_vec(U)
#
# mat"
# figure
# hold on
# plot($um')
# set(gca, 'YScale', 'log')
# hold off
# "

U1 = deepcopy(U)
A,B = getAB(model,X,U,dt)
@assert (length(A) == length(B))
@assert (length(X) == length(U)+1)
@assert (length(A) == length(U))
N = length(X)
x0 = copy(X[1])

# jmodel = Model(Mosek.Optimizer)
# set_optimizer_attribute(jmodel, "MSK_DPAR_INTPNT_CO_TOL_DFEAS",1e-15)
# set_optimizer_attribute(jmodel, "MSK_DPAR_INTPNT_CO_TOL_REL_GAP",1e-15)

jmodel = Model(COSMO.Optimizer)
# set_optimizer_attribute(jmodel, "eps_abs",1e-8)


@variable(jmodel, δx[1:6,1:N])
@variable(jmodel, δu[1:2,1:N-1])

@constraint(jmodel,δx[:,1] .== zeros(6))
for i = 1:N-1
    @constraint(jmodel, δx[:,i+1] .== A[i]*δx[:,i] + B[i]*δu[:,i])
end

# get maximum lift
maxL = zeros(N-1)
for i = 1:N-1
    maxL[i] = getmaxL(model,X[i])
end

for i = 1:N-1
    @constraint(jmodel, [maxL[i],δu[1,i],δu[2,i]] in SecondOrderCone())
end

rr = normalize(xf[1:3])
Qn = I - rr*rr'
Q = Qn'*Qn
Q = I(3)
objective_exp = @expression(jmodel, (X[N][1:3] + δx[1:3,N] - xf[1:3] )' * Q * (X[N][1:3] + δx[1:3,N] - xf[1:3] ))

α = 1e-2
for i = 1:N-1
    add_to_expression!(objective_exp, α*δu[:,i]'*δu[:,i])
end

@objective(jmodel, Min, objective_exp)

optimize!(jmodel)

du = value.(δu)
dx = value.(δx)
mat"
figure
hold on
plot($du')
hold off
"

U = vec_from_mat(mat_from_vec(U) + du)

dt = 2/3600
N = 160
X = NaN*[@SArray zeros(6) for i = 1:N]
# U = [@SArray zeros(2) for i = 1:N-1]

X[1] = deepcopy(x0)
# global end_idx = NaN
for i = 1:(N-1)
    # U[i] = getmaxL(model,X[i])*[0;.5]
    X[i+1] = rk4(model,X[i],U[i],dt)
    if altitude(model,X[i+1])<10
        @info "under altitude"
        # end_idx = i +1
        # X = X[1:(i+1)];U = U[1:(i+1)];
        break
    end
end

# xf_dr, xf_cr = rangedistances(model,xf,x0)
alt, dr, cr = postprocess(model::EntryVehicle,X,x0)
mat"
figure
hold on
plot($dr,$cr)
plot($xf_dr,$xf_cr,'r.','markersize',20)
xlabel('Downrange')
ylabel('Crossrange')
legend('Trajectory','Goal','location','northwest')
hold off
"

u1m = mat_from_vec(U1)
um = mat_from_vec(U)
mat"
figure
hold on
plot($u1m','b')
plot($um','r')
hold off
"

mat"
figure
hold on
plot($alt)
plot(1:length($alt),ones( length($alt),1)*10,'r' )
hold off
"
