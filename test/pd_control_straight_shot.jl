const RD = RobotDynamics
using StaticArrays
using TrajectoryOptimization
const TO = TrajectoryOptimization
using EntryGuidance
const EG = EntryGuidance
using Random


struct EntryVehicle_fixed_time_direct{T} <: TO.AbstractModel
    evmodel::EG.CartesianModel{T}
end
function RD.dynamics(model::EntryVehicle_fixed_time_direct, x, u)
    # U here is [D;L;σ] where σ is bank angle
    return EG.dynamics(x[1:6], u, model.evmodel)
end

# Define a custom integration method
abstract type EntryVehicleRK_fixed_time_direct <: RD.Explicit end

# Define the discrete dynamics function
function RD.discrete_dynamics(::Type{EntryVehicleRK_fixed_time_direct}, model::EntryVehicle_fixed_time_direct,
        x::StaticVector, u::StaticVector, t, dt)
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

N = 2000
dt = 0.1
tf = 250
t_vec = 0:dt:tf
N = length(t_vec)
X = NaN*[@SVector zeros(6) for i = 1:N]
U = NaN*[@SVector zeros(2) for i = 1:N]
tscale = 3600
x0 = [r0;v0]
Random.seed!(123)
X[1] = x0 #+ SVector{6}([1*normalize(randn(3));.0001*tscale*normalize(randn(3))])
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
    # this sim U is [D,L,σ]
    α = deg2rad(15)
    β = 0
    sim_u = EG.angles_input([α, β],X[i][1:6],model.evmodel)
    z = KnotPoint(X[i],sim_u,dt)
    X[i+1] = discrete_dynamics(EntryVehicleRK_fixed_time_direct,model,z)
end

# Now we do this a second time
X2 = NaN*[@SVector zeros(6) for i = 1:N]
U2 = NaN*[@SVector zeros(2) for i = 1:N]
# Random.seed!(123)
X2[1] = x0 + SVector{6}([1*normalize(randn(3));.0001*tscale*normalize(randn(3))])


α_hist = zeros(2000)
β_hist = zeros(2000)

for i = 1:2000#(N-1000)

    # find closest point
    ran = 1:3
    dist = ([norm(X2[i][ran] - X[k][ran]) for k = 1:length(X)])
    for i = 1:length(dist)
        if isnan(dist[i])
            dist[i] = Inf
        end
    end
    idx = argmin(dist)
    # get position errors
    re_hist[i], ve_hist[i] = get_lvlh_errors(X2[i][1:3],X2[i][4:6],X[idx][1:3],X[idx][4:6])

    # this sim U is [D,L,σ]
    kp = 50
    kd = 0#.05
    α = clamp(deg2rad(10) -kp*re_hist[i][2] - kd*ve_hist[i][2],-deg2rad(25),deg2rad(25))
    # α = deg2rad(10)
    kp = 50
    kd = 0#.15
    β = clamp(deg2rad(-kp*re_hist[i][1] -kd*ve_hist[i][1]),-deg2rad(60),deg2rad(60))

    α_hist[i] = α
    β_hist[i] = β
    sim_u = EG.angles_input([α, β],X2[i][1:6],model.evmodel)
    z = KnotPoint(X2[i],sim_u,dt)
    # @infiltrate
    # error()
    X2[i+1] = discrete_dynamics(EntryVehicleRK_fixed_time_direct,model,z)
end

using MATLAB
using Attitude
xm_altro = mat_from_vec(X2)
xm_sim = mat_from_vec(X)
# t_traj = traj.t
# mat"
#
# figure
# hold on
# sgtitle('MCMF Position')
# subplot(3,1,1)
# hold on
# plot($t_vec, $xm_altro(1,:))
# plot($t_vec,$xm_sim(1,:))
# legend('Controlled','Uncontrolled')
# subplot(3,1,2)
# hold on
# plot($t_vec, $xm_altro(2,:))
# plot($t_vec,$xm_sim(2,:))
# legend('Controlled','Uncontrolled')
# subplot(3,1,3)
# hold on
# plot($t_vec, $xm_altro(3,:))
# plot($t_vec,$xm_sim(3,:))
# legend('Controlled','Uncontrolled')
# hold off
# "
# mat"
# figure
# sgtitle('MCMF Velocity')
# hold on
# subplot(3,1,1)
# hold on
# plot($t_vec, $xm_altro(4,:))
# plot($t_vec,$xm_sim(4,:))
# legend('Controlled','Uncontrolled')
# subplot(3,1,2)
# hold on
# plot($t_vec, $xm_altro(5,:))
# plot($t_vec,$xm_sim(5,:))
# legend('Controlled','Uncontrolled')
# subplot(3,1,3)
# hold on
# plot($t_vec, $xm_altro(6,:))
# plot($t_vec,$xm_sim(6,:))
# legend('Controlled','Uncontrolled')
# hold off
# "

rem = mat_from_vec(re_hist)
vem = mat_from_vec(ve_hist)
mat"
figure
hold on
plot($rem(1:2,:)')
legend('X','Y')
hold off
"
mat"
figure
hold on
subplot(2,1,1)
plot(rad2deg($α_hist))
legend('Alpha')

subplot(2,1,2)
plot(rad2deg($β_hist))
legend('Beta')
hold off
"
# mat"
# figure
# hold on
# plot($rem(1,:),$rem(2,:))
# hold off
# grid on
# axis equal
# "

#
#
# um_altro = mat_from_vec(Up)
# um_sim = mat_from_vec(U)
#
# mat"
# figure
# hold on
# subplot(2,1,1)
# hold
# title('Bank Angle Derivative')
# plot($t_traj(1:end),$um_altro(1,:))
# plot($t_vec,$um_sim(1,:))
#
# subplot(2,1,2)
# title('Angle of Attack')
# hold on
# plot($t_traj(1:end),$um_altro(2,:))
# plot($t_vec,$um_sim(2,:))
#
# hold off
# "
# rem = mat_from_vec(re_hist)
# vem = mat_from_vec(ve_hist)
# mat"
# figure
# hold on
# plot($rem(1:2,:)')
# hold off
# "
# mat"
# figure
# hold on
# plot($rem(1,:),$rem(2,:))
# hold off
# grid on
# "
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
