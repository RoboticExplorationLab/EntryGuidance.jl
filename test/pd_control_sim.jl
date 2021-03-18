const RD = RobotDynamics
using StaticArrays

struct EntryVehicle_fixed_time_direct{T} <: TO.AbstractModel
    evmodel::EG.CartesianModel{T}
end
function RD.dynamics(model::EntryVehicle_fixed_time_direct, x, u)
    # α = u[2]
    # σ̇ = u[1]
    # return [EG.dynamics(x[1:6], EG.angles_input([α, x[7]],x[1:6],model.evmodel), model.evmodel); σ̇]
    # u = EG.angles_input([α, x[7]],x[1:6],model.evmodel)
    return EG.dynamics(x[1:6], u, model.evmodel)
end

# Define a custom integration method
abstract type EntryVehicleRK_fixed_time_direct <: RD.Explicit end

# Define the discrete dynamics function
function RD.discrete_dynamics(::Type{EntryVehicleRK_fixed_time_direct}, model::EntryVehicle_fixed_time_direct,
        x::StaticVector, u::StaticVector, t, dt)

    # h = u[3]/3600.0 #u is in seconds, dynamics are in hours
    h = dt/3600 # dt is seconts, h is in hours

    k1 = RD.dynamics(model, x,             u)*h;
    # @infiltrate
    # error()
    k2 = RD.dynamics(model, x + k1/2,      u)*h;
    k3 = RD.dynamics(model, x - k1 + 2*k2, u)*h;

    return x + (k1 + 4*k2 + k3)/6
end

Base.size(::EntryVehicle_fixed_time_direct) = 6,2
model = EntryVehicle_fixed_time_direct(CartesianMSLModel())
n,m = size(model)

Xp = traj.X
Up = traj.U
N = length(X)
dt = 0.1
tf = 250
t_vec = 0:dt:tf
N = length(t_vec)
X = NaN*[@SVector zeros(6) for i = 1:N]
U = NaN*[@SVector zeros(2) for i = 1:N]
tscale = 3600
x0 = [r0;v0]
Random.seed!(123)
X[1] = x0 + SVector{6}([1*normalize(randn(3));.0001*tscale*normalize(randn(3))])
re_hist = NaN*[zeros(3) for i = 1:N]
ve_hist = NaN*[zeros(3) for i = 1:N]
function get_lvlh_errors(r,v,rx,vx)
    e1 = normalize(cross(r,v))
    e2 = normalize(cross(v,e1))
    e3 = normalize(cross(e1,e2))

    ECI_Q_LVLH = [e1 e2 e3]

    # position error
    r_e = (ECI_Q_LVLH')*(r - rx)
    v_e = (ECI_Q_LVLH')*(v - vx)

    return r_e, v_e
end
for i = 1:2000#(N-1000)
    x = copy(X[i])
    # find closest point on trajectory
    ran = 1:3
    dist = [norm(x[ran] - Xp[k][ran]) for k = 1:length(X)]
    idx = argmin(dist)

    # errors
    # @infiltrate
    # error()
    re_hist[i], ve_hist[i] = get_lvlh_errors(Xp[idx][1:3],Xp[idx][4:6],X[i][1:3],X[i][4:6])
    # @infiltrate
    # error()
    # control plan
    # U[i] = Up[idx][1:2]

    # step forward in simulation
    # this sim U is [D,L,σ]
    α = 1*Up[idx][2]
    # if i < 1000 && i > 800
        kp = 1
        kd = 1
        δβ = atan(re_hist[i][1],re_hist[i][2]) + atan(ve_hist[i][1],ve_hist[i][2])
        # β = Xp[idx][7]+δβ
        β = copy(δβ)
        @show rad2deg(δβ)
    # else
    #     β = Xp[idx][7]
    # end
    # β = Xp[idx][7]
    U[i] = [α;β]
    # @show rad2deg(α)
    # if mod(i,50) == 0
    #     ρ = atmospheric_density(X[i][1:3], model.evmodel)
    #     # @show ρ
    # end
    sim_u = EG.angles_input([α, β],X[i][1:6],model.evmodel)
    sim_u[3]   = β
    # @show rad2deg(α)
    # @show rad2deg(sim_u[3])
    # if i < 1000 && i > 800
        # sim_u[3] = β+atan(re_hist[i][1],re_hist[i][2])
    # else
    #     sim_u[3] = β
    # end
    # @infiltrate
    # error()
    z = KnotPoint(X[i],sim_u,dt)
    # @infiltrate
    # error()
    X[i+1] = discrete_dynamics(EntryVehicleRK_fixed_time_direct,model,z)
end



using MATLAB
using Attitude
xm_altro = mat_from_vec(Xp)
xm_sim = mat_from_vec(X)
t_traj = traj.t
mat"

figure
hold on
sgtitle('MCMF Position')
subplot(3,1,1)
hold on
plot($t_traj, $xm_altro(1,:))
plot($t_vec,$xm_sim(1,:))
subplot(3,1,2)
hold on
plot($t_traj, $xm_altro(2,:))
plot($t_vec,$xm_sim(2,:))
subplot(3,1,3)
hold on
plot($t_traj, $xm_altro(3,:))
plot($t_vec,$xm_sim(3,:))
hold off
"
mat"
figure
sgtitle('MCMF Velocity')
hold on
subplot(3,1,1)
hold on
plot($t_traj, $xm_altro(4,:))
plot($t_vec,$xm_sim(4,:))
subplot(3,1,2)
hold on
plot($t_traj, $xm_altro(5,:))
plot($t_vec,$xm_sim(5,:))
subplot(3,1,3)
hold on
plot($t_traj, $xm_altro(6,:))
plot($t_vec,$xm_sim(6,:))
hold off
"
#
#
um_altro = mat_from_vec(Up)
um_sim = mat_from_vec(U)
#
mat"
figure
hold on
subplot(2,1,1)
hold
title('Bank Angle Derivative')
%plot($t_traj(1:end),$um_altro(1,:))
plot($t_vec,$um_sim(1,:))

subplot(2,1,2)
title('Angle of Attack')
hold on
%plot($t_traj(1:end),$um_altro(2,:))
plot($t_vec,$um_sim(2,:))

hold off
"
rem = mat_from_vec(re_hist)
vem = mat_from_vec(ve_hist)
mat"
figure
hold on
plot($rem(1:2,:)')
hold off
"
mat"
figure
hold on
plot($rem(1,:),$rem(2,:))
hold off
grid on
axis equal
"
# difrem = diff(rem;dims = 2)/dt
# mat"
# figure
# hold on
# plot($difrem(1:2,:)')
# hold off
# "
# mat"
# figure
# hold on
# plot($vem(1:2,:)')
# hold off
# "


Xpm = mat_from_vec(Xp)
mat"
figure
hold on
plot($Xpm(7,:))
hold off
"
