using LinearAlgebra, Attitude, JLD2, MATLAB, RobotDynamics
using TrajectoryOptimization
const TO = TrajectoryOptimization
using EntryGuidance
const EG = EntryGuidance

@load "goodtraj.jld2" traj

struct EntryVehicle{T} <: TO.AbstractModel
    evmodel::EG.CartesianModel{T}
end
model = EntryVehicle(CartesianMSLModel())

# initial conditions
Rm = model.evmodel.planet.R
r0 = [Rm+125.0, 0.0, 0.0] #Atmospheric interface at 125 km altitude
V0 = 5.845*3600 #Mars-relative velocity at interface of 5.845 km/sec
γ0 = -15.474*(pi/180.0) #Flight path angle at interface
v0 = V0*[sin(γ0), cos(γ0), 0.0]

function post_process_traj(traj)
    X = traj.X
    U = traj.U
    t = traj.t
    dt = traj.dt
    Rm = 3396.2

    N = length(X)
    alt = zeros(N)
    bank = zeros(N)
    #
    for k = 1:N
        alt[k] = norm(X[k][1:3])-Rm
        bank[k] = X[k][7]#x_traj[7,k]
    end

    # u_traj = zeros(m,N-1)
    lift = zeros(2,N-1)
    AoA = zeros(N-1)
    σ̇ = zeros(N-1)
    down_range = zeros(N)
    cross_range = zeros(N)
    # @infiltrate
    # error()
    for k = 1:(N-1)

        # get the controls stuff
        D,L,σ = EG.angles_input([U[k][2], X[k][7]],X[k][1:6],model.evmodel)
        r = X[k][1:3]
        v = X[k][4:6]
        e1 = cross(r,v)
        e1 = e1./norm(e1)
        e2 = cross(v,e1)
        e2 = e2./norm(e2)
        a_d = -(D/norm(v)).*v
        # a_l = L*sin(σ)*e1 + L*cos(σ)*e2
        lift[:,k] = [L*sin(σ); L*cos(σ)]

        σ̇[k] = U[k][1]#u_traj[1,k]
        AoA[k] = U[k][2]#u_traj[2,k]
        down_range[k+1] = down_range[k] + norm(X[k+1][1:3] - X[k][1:3])
        cross_range[k+1] = cross_range[k] + (cross(r0,v0)/norm(cross(r0,v0)))'*(X[k+1][1:3] - X[k][1:3])
    end

    # bank_angle = copy(bank)
    return lift, bank, AoA

end

lift2, bank2, AoA2 = post_process_traj(traj)

mat"
figure
hold on
plot($lift2(1,:),$lift2(2,:))
plot($lift2(1,1),$lift2(2,1),'ro')
plot($lift2(1,end),$lift2(2,end),'go')
legend('Lift Vector','Start','Stop')
hold off
"
mat"
figure
subplot(2,1,1)
plot($lift2(1,:))

subplot(2,1,2)
plot($lift2(2,:))
"
mat"
figure
hold on
title('Bank Angle')
plot($bank2)
hold off
"

mat"
figure
hold on
title('Angle of Attack')
plot($AoA2)
hold off
"


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


get_lvlh_errors(r0,v0,0,0)
