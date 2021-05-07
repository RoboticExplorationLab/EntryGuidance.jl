using EntryGuidance
const EG = EntryGuidance
using ForwardDiff
using LinearAlgebra
using MATLAB
using StaticArrays
using Attitude
# evmod = CartesianMSLModel()
struct EntryVehicle{T}
    evmodel::EG.CartesianModel{T}
end
model = EntryVehicle(CartesianMSLModel())
function evdynamics(M::EntryVehicle, x, u)
    α = u[1]
    σ = u[2]
    return EG.dynamics(x[1:6], EG.angles_input([α, σ],x[1:6],M.evmodel), M.evmodel)
end
function rk4(model::EntryVehicle,x_n,u,dt)
    k1 = dt*evdynamics(model,x_n,u)
    k2 = dt*evdynamics(model,x_n+k1/2,u)
    k3 = dt*evdynamics(model,x_n+k2/2,u)
    k4 = dt*evdynamics(model,x_n+k3,u)
    return (x_n + (1/6)*(k1 + 2*k2 + 2*k3 + k4))
end
function plot(x::Vector,y::Vector)
    mat"
    figure
    hold on
    plot($x,$y)
    hold off
    "
end
function plot(y::Vector,title::String)
    mat"
    figure
    hold on
    title($title)
    plot($y)
    hold off
    "
end
function plot(y)
    mat"
    figure
    hold on
    plot($y)
    hold off
    "
end
function get_heating(model::EntryVehicle, x)

    r = x[1:3]

    # # convert v to m/s
    # v = x[4:6]*1000/3600
    # # v = x[4:6]*1/3600
    #
    # #atmospheric density # units of kg / m^3
    # ρ = atmospheric_density(r, model.evmodel)*1e-9
    #
    # # heatrate W/m^2
    # HR = 8.43e-13*ρ^(0.82958)*norm(v)^4.512

    # convert v to m/s
    v = x[4:6]*1000/3600

    #atmospheric density # units of kg / m^3
    ρ = atmospheric_density(r, model.evmodel)*1e-9
    # ρ = myρ(x)
    # ρ*= 1e6

    # heatrate W/m^2
    HR = 8.43e-13*ρ^(0.82958)*norm(v)^4.512
    # rn = 0.6 # meters
    # rn = 1
    # κ = 1.9027e-4
    #
    # HR = κ*sqrt(ρ/rn)*norm(v)^3 # per m^3
    # HR *= 1e-6

    # pressure
    # PR = 0.80527*ρ^(1.0036)*norm(v)^2.0251
    PR = get_pressure(model,x)

    return HR, PR
end
function get_pressure(model::EntryVehicle,x)
    # input is the km, km/s units
    # output is the scaled kPa units
    ρ = atmospheric_density(x[1:3], model.evmodel)*1e-9
    # PR = 0.80527*ρ^(1.0036)*norm(x[4:6]*1000/3600)^2.0251
    PR = 0.5*ρ*norm(x[4:6]*1000/3600)^2
    PR /= 1000
    return PR
end
function safety_post(model::EntryVehicle,X,dt)

    N = length(X)
    HR_hist = zeros(N)
    PR_hist = zeros(N)
    g_hist = zeros(N-1)

    for i = 1:N
        HR_hist[i],PR_hist[i] = get_heating(model::EntryVehicle, X[i])
    end
    for i = 1:N-1
        v0 = X[i][4:6]*1000/3600
        v1 = X[i+1][4:6]*1000/3600
        dt = 2# seconds
        g_hist[i] = norm(v1 - v0)/(dt*9.8)
    end
    return HR_hist, PR_hist, g_hist
end

Rm = model.evmodel.planet.R
r0 = [Rm+125.0, 0.0, 0.0] #Atmospheric interface at 125 km altitude
V0 = 5.845*3600 #Mars-relative velocity at interface of 5.845 km/sec
γ0 = -15.474*(pi/180.0) #Flight path angle at interface
v0 = V0*[sin(γ0), cos(γ0), 0.0]
α = 15.0*pi/180 #derived from MSL initial L/D = 0.24
σ = 0.0*pi/180 #this is totally made up
x0 = [r0; v0]
function rollout(model,x0,U_in,dt)
    N = 1000
    X = [zeros(6) for i = 1:N]
    U = [zeros(2) for i = 1:N-1]
    X[1] = copy(x0)
    end_idx = NaN
    t_vec = 0:dt:((N-1)*dt)
    for i = 1:(length(X)-1)

        U[i] = i > length(U_in) ? deepcopy(U_in[end]) : deepcopy(U_in[i])

        # step forward in sim
        X[i+1] = rk4(model,X[i],U[i],dt)

        # check for hitting the 10,000ft mark
        if (norm(X[i+1][1:3]) -  model.evmodel.planet.R) < 10
            end_idx = i+1
            @info "rollout successful"
            break
        end
    end
    # trim relavent
    if isnan(end_idx)
        error("didn't hit the altitude during the rollout")
    end

    X = X[1:end_idx]
    U = U[1:(end_idx-1)]
    return X, U
end

function firsttest()
U = [[deg2rad(18),0] for i = 1:100]
dt = 2/3600.0
X2, U2 = rollout(model, x0, U, dt )

HR_hist, PR_hist, g_hist = safety_post(model::EntryVehicle,X2,dt)

# plot(HR_hist, "Heat Rate")
# plot(PR_hist, "Pressure kPa")
# plot(g_hist, "G hist")


H = [zeros(6) for i = 1:length(X2)]
for i = 1:length(H)
    hess = ForwardDiff.hessian(_x -> get_pressure(model::EntryVehicle,_x),X2[i])
    H[i] = eigvals(hess)
end

H = mat_from_vec(H)
mat"
figure
hold on
plot($H')
hold off
"
end

firsttest()
