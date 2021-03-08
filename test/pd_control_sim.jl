const RD = RobotDynamics
using StaticArrays

struct EntryVehicle_fixed_time{T} <: TO.AbstractModel
    evmodel::EG.CartesianModel{T}
end
function RD.dynamics(model::EntryVehicle_fixed_time, x, u)
    α = u[2]
    σ̇ = u[1]
    return [EG.dynamics(x[1:6], EG.angles_input([α, x[7]],x[1:6],model.evmodel), model.evmodel); σ̇]
end

# Define a custom integration method
abstract type EntryVehicleRK_fixed_time <: RD.Explicit end

# Define the discrete dynamics function
function RD.discrete_dynamics(::Type{EntryVehicleRK_fixed_time}, model::EntryVehicle_fixed_time,
        x::StaticVector, u::StaticVector, t, dt)

    # h = u[3]/3600.0 #u is in seconds, dynamics are in hours
    h = dt/3600 # dt is seconts, h is in hours

    k1 = RD.dynamics(model, x,             u)*h;
    k2 = RD.dynamics(model, x + k1/2,      u)*h;
    k3 = RD.dynamics(model, x - k1 + 2*k2, u)*h;

    return x + (k1 + 4*k2 + k3)/6
end

Base.size(::EntryVehicle_fixed_time) = 7,2
model = EntryVehicle_fixed_time(CartesianMSLModel())
n,m = size(model)

Xp = traj.X
Up = traj.U
N = length(X)
dt = 0.1
tf = 250
t_vec = 0:dt:tf
N = length(t_vec)
X = NaN*[@SVector zeros(7) for i = 1:N]
U = NaN*[@SVector zeros(2) for i = 1:N]
tscale = 3600
x0 = [r0;v0;deg2rad(15)]
X[1] = x0 + SVector{7}([10*normalize(randn(3));.0001*tscale*normalize(randn(3));0])
# X[1] = x0
re_hist = NaN*[zeros(3) for i = 1:N]
# function get_lvlh_errors(r,v,rx,vx)
#     e1 = normalize(cross(r,v))
#     e2 = normalize(cross(v,e1))
#     e3 = normalize(cross(e1,e2))
#
#     ECI_Q_LVLH = [e1 e2 e3]
#
#     # position error
#     r_e = (ECI_Q_LVLH')*(r - rx)
#     v_e = (ECI_Q_LVLH')*(v - vx)
#
#     return r_e, v_e
# end
for i = 1:2000#(N-1000)
    x = copy(X[i])
    # find closest point on trajectory
    ran = 1:3
    dist = [norm(x[ran] - Xp[k][ran]) for k = 1:length(X)]
    idx = argmin(dist)

    # errors
    # @infiltrate
    # error()
    re_hist[i], v_e = get_lvlh_errors(Xp[idx][1:3],Xp[idx][4:6],X[i][1:3],X[i][4:6])
    # @infiltrate
    # error()
    # control plan
    U[i] = Up[idx][1:2]

    # step forward in simulation
    z = KnotPoint(X[i],U[i],dt)
    X[i+1] = discrete_dynamics(EntryVehicleRK_fixed_time,model,z)
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


um_altro = mat_from_vec(Up)
um_sim = mat_from_vec(U)

mat"
figure
hold on
subplot(2,1,1)
hold
title('Bank Angle Derivative')
plot($t_traj(1:end),$um_altro(1,:))
plot($t_vec,$um_sim(1,:))

subplot(2,1,2)
title('Angle of Attack')
hold on
plot($t_traj(1:end),$um_altro(2,:))
plot($t_vec,$um_sim(2,:))

hold off
"
rem = mat_from_vec(re_hist)

mat"
figure
hold on
plot($rem(1:2,:)')
hold off
"
