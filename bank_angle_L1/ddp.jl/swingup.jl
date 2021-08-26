using ForwardDiff, StaticArrays, MATLAB, Infiltrator, LinearAlgebra
# struct SAT
#     J::Diagonal{Float64,SArray{Tuple{3},Float64,1,3}}
# end
#
# struct LQR
#     Q::Diagonal{Float64,SArray{Tuple{6},Float64,1,6}}
#     R::Diagonal{Float64,SArray{Tuple{3},Float64,1,3}}
#     Qf::Diagonal{Float64,SArray{Tuple{6},Float64,1,6}}
#     xf::SArray{Tuple{6},Float64,1,6}
# end
# # fill in parameters
# # EOM() = EOM(4,1,
# #     1.0, 0.2, 0.5, 9.81)
# model = SAT(Diagonal(SA[1,2,3.0]))
# # define equations of motion
# function dynamics(model::SAT, x, u)
#     J = model.J
#     p = x[SA[1,2,3]]
#     ω = x[SA[4,5,6]]
#     α = J\(u - cross(ω,J*ω))
#     pdot = pdot_from_w(p,ω)
#     return SA[pdot[1],pdot[2],pdot[3],α[1],α[2],α[3]]
# end
# function rk4(model,x_n,u,h)
#     k1 = h*dynamics(model,x_n,u)
#     k2 = h*dynamics(model,x_n+k1/2,u)
#     k3 = h*dynamics(model,x_n+k2/2,u)
#     k4 = h*dynamics(model,x_n+k3,u)
#     return (x_n + (1/6)*(k1+2*k2+2*k3 + k4))
# end
#
# function getAB(model,x,u,dt)
#     A = ForwardDiff.jacobian(_x -> rk4(model,_x,u,dt),x)
#     B = ForwardDiff.jacobian(_u -> rk4(model,x,_u,dt),u)
#     return A, B
# end
#
# function cost_expansion(cost::LQR,x,u)
#     Q = deepcopy(cost.Q)
#     q = cost.Q*(x - cost.xf)
#     R = deepcopy(cost.R)
#     r = cost.R*u
#     return Q, q, R, r
# end
# function evalcost(cost::LQR,X,U)
#     N = length(X)
#     J = 0.0
#     Q, R, xf, Qf = cost.Q, cost.R, cost.xf, cost.Qf
#     for i = 1:length(X-1)
#         x = X[i]
#         u = U[i]
#         J += (x - xf)'*Q*(x - xf) + u'*R*u
#     end
#     J += (X[N]- xf)'*Qf*(X[N]-xf)
#     return J
# end
# function naiverollout(model,X,U,dt,x0)
#     X[1] = deepcopy(x0)
#     N = length(X)
#     for i = 1:(N-1)
#         X[i+1] = rk4(model,X[i],U[i],dt)
#     end
#     return X
# end

# function backwards_pass(model,X,U,dt)
#
#
#
function Tmat(n)
    n2 = n^2
    T = spzeros(n2,n2)
    for i = 1 : n
        for j = 1 : n
            T[i + n*(j - 1), j + n*(i - 1)] = 1
        end
    end
    return T
end
function Avec(x,u)
    return vec(ForwardDiff.jacobian(dx->dynamics_rk4(dx,u),x))
end
function Bvec(x,u)
    return vec(ForwardDiff.jacobian(du->dynamics_rk4(x,du),u))
end
const Jinertia = Diagonal([1,2,3])
function dynamics(x, u)
    J = Jinertia
    p = x[SA[1,2,3]]
    ω = x[SA[4,5,6]]
    # @infiltrate
    α = J\(u - cross(ω,J*ω))
    pdot = pdot_from_w(p,ω)
    return [pdot[1],pdot[2],pdot[3],α[1],α[2],α[3]]
end
function dynamics_rk4(x,u)
    #RK4 integration with zero-order hold on u
    f1 = dynamics(x, u)
    f2 = dynamics(x + 0.5*h*f1, u)
    f3 = dynamics(x + 0.5*h*f2, u)
    f4 = dynamics(x + h*f3, u)
    return x + (h/6.0)*(f1 + 2*f2 + 2*f3 + f4)
end
h = .2
dt = h
Nx = 6     # number of state
Nu = 3     # number of controls
Tfinal = 20.0 # final time
Nt = Int(Tfinal/h)+1    # number of time steps
thist = Array(range(0,h*(Nt-1), step=h));
# Cost weights
Q = Diagonal(ones(Nx))
R = Diagonal(ones(Nu))
Qn = 10*Q
# Q = Diagonal([1.0*ones(2); 0.1*ones(2)]);
# R = 0.5;
# Qn = Array(1000.0*I(Nx));
function stage_cost(x,u)
    return 0.5*((x-xgoal)'*Q*(x-xgoal)) + 0.5*u'*R*u
end
function terminal_cost(x)
    return 0.5*(x-xgoal)'*Qn*(x-xgoal)
end
function cost(xtraj,utraj)
    J = 0.0
    for k = 1:(Nt-1)
        J += stage_cost(xtraj[k],utraj[k])
    end
    J += terminal_cost(xtraj[Nt])
    return J
end
#Initial guess
# x0 = [-pi/2; 0; 0; 0]


function runme(ddp)
    p = [zeros(Nx) for i = 1:Nt]
# P = zeros(Nx,Nx,Nt)
P = [zeros(Nx,Nx) for i = 1:Nt]
d = [zeros(Nu) for i = 1:(Nt-1)]
# K = zeros(Nu,Nx,Nt-1)
K = [zeros(Nu,Nx) for i = 1:Nt-1]
ΔJ = 0.0

x0 = [p_from_phi(deg2rad(170)*normalize([1;2;3])); deg2rad(0)*normalize([1;2;3])]
xgoal = zeros(6)
# xtraj = kron(ones(1,Nt), x0)
xtraj = [zeros(length(x0)) for i = 1:Nt]
xtraj[1] = x0
# utraj = .1*randn(Nt-1);
utraj = [.1*zeros(Nu) for i = 1:(Nt-1)]
#Initial Rollout
for k = 1:(Nt-1)
    xtraj[k+1] .= dynamics_rk4(xtraj[k],utraj[k])
end
J = cost(xtraj,utraj)

gx = zeros(Nx)
gu = zeros(Nu)
Gxx = zeros(Nx,Nx)
Guu = zeros(Nu,Nu)
Gxu = zeros(Nx,Nu)
Gux = zeros(Nu,Nx)

iter = 0
for i = 1:100
    iter += 1
    ΔJ = 0.0


    p[Nt] = Qn*(xtraj[Nt]-xgoal)
    P[Nt] = Qn

    #Backward Pass
    for k = (Nt-1):-1:1
        #Calculate derivatives
        q = Q*(xtraj[k]-xgoal)
        r = R*utraj[k]

        A = ForwardDiff.jacobian(dx->dynamics_rk4(dx,utraj[k]),xtraj[k])
        B = ForwardDiff.jacobian(du->dynamics_rk4(xtraj[k],du),utraj[k])

        gx = q + A'*p[k+1]
        gu = r + B'*p[k+1]

        if ddp
            dAdx = ForwardDiff.jacobian(_x -> Avec(_x,utraj[k]),xtraj[k])
            dBdx = ForwardDiff.jacobian(_x -> Bvec(_x,utraj[k]),xtraj[k])
            dBdu = ForwardDiff.jacobian(_u -> Bvec(xtraj[k],_u),utraj[k])

            fxx = reshape(dAdx,6,6,6)
            fuu = reshape(dBdu,6,3,3)
            fux = reshape(dBdx,6,3,6)

            Gxx .= Q + A'*P[k+1]*A + p[k+1]' ⊡ fxx
            Guu .= R + B'*P[k+1]*B + p[k+1]' ⊡ fuu
            Gux .= B'*P[k+1]*A     + p[k+1]' ⊡ fux
            Gxu .= copy(transpose(Gux))
        else
            Gxx .= Q + A'*P[k+1]*A #+ p[k+1]' ⊡ fxx
            Guu .= R + B'*P[k+1]*B #+ p[k+1]' ⊡ fuu
            Gux .= B'*P[k+1]*A #+ p[k+1]' ⊡ fux
            Gxu .= copy(transpose(Gux))
        end

        ρ = 1e-5
        d[k] = (Guu+ρ*I)\gu
        K[k] = (Guu+ρ*I)\Gux

        p[k] = gx - K[k]'*gu + K[k]'*Guu*d[k] - Gxu*d[k]
        P[k] = Gxx + K[k]'*Guu*K[k] - Gxu*K[k] - K[k]'*Gux

        ΔJ += gu'*d[k]
    end


    #Forward rollout with line search
    xn = copy(xtraj)
    un = copy(utraj)
    xn[1] = copy(xtraj[1])
    α = 1.0
    Jn = Inf

    for kk = 1:20
        xn = copy(xtraj)
        un = copy(utraj)
        for k = 1:(Nt-1)
            un[k] = utraj[k] - α*d[k] - K[k]*(xn[k]-xtraj[k])
            xn[k+1] = dynamics_rk4(xn[k],un[k])
        end
        Jn = cost(xn,un)
        if Jn < (J)# - 1e-2*α*ΔJ)
            @info "α = $α"
            break
        else
            α *= .5
        end
        if kk == 20
            error("line search failed")
        end

    end
    @show iter
    dJ = abs(J - Jn)
    J = copy(Jn)
    @show dJ
    xtraj .= xn
    utraj .= un
    if dJ<1e-8
        break
    end
end

xm = mat_from_vec(xtraj)
mat"
figure
hold on
plot($xm')
hold off
"
end

runme(false)
