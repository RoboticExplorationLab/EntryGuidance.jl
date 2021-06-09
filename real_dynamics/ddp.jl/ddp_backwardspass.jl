using TensorCore

function Avec(model,x,u,dt)
    return vec(ForwardDiff.jacobian(dx->rk4(model,dx,u,dt),x))
end
function Bvec(model,x,u,dt)
    return vec(ForwardDiff.jacobian(du->rk4(model,x,du,dt),u))
end

function backpass(model,xtraj,utraj,xgoal,dt,ddp)
Nt = length(xtraj)
Nx = length(xtraj[1])
Nu = length(utraj[1])
p = [zeros(Nx) for i = 1:Nt]
P = [zeros(Nx,Nx) for i = 1:Nt]
d = [zeros(Nu) for i = 1:(Nt-1)]
K = [zeros(Nu,Nx) for i = 1:Nt-1]

Q = 1*Diagonal(@SVector ones(6))
Qn = Diagonal(@SVector ones(6))
R = 1*Diagonal(@SVector ones(2))

p[Nt] = zeros(Nx)#Qn*(xtraj[Nt]-xgoal)
P[Nt] = Qn

#Backward Pass
for k = (Nt-1):-1:1
    #Calculate derivatives
    # q = Q*(xtraj[k]-xgoal)
    # r = R*utraj[k]
    q = zeros(Nx)
    r = zeros(Nu)

    A = ForwardDiff.jacobian(dx->rk4(model,dx,utraj[k],dt),xtraj[k])
    B = ForwardDiff.jacobian(du->rk4(model,xtraj[k],du,dt),utraj[k])

    gx = q + A'*p[k+1]
    gu = r + B'*p[k+1]

    if ddp
        dAdx = ForwardDiff.jacobian(_x -> Avec(model,_x,utraj[k],dt),xtraj[k])
        dBdx = ForwardDiff.jacobian(_x -> Bvec(model,_x,utraj[k],dt),xtraj[k])
        dBdu = ForwardDiff.jacobian(_u -> Bvec(model,xtraj[k],_u,dt),utraj[k])

        fxx = reshape(dAdx,6,6,6)
        fuu = reshape(dBdu,6,2,2)
        fux = reshape(dBdx,6,2,6)

        Gxx = Q + A'*P[k+1]*A + p[k+1]' ⊡ fxx
        Guu = R + B'*P[k+1]*B + p[k+1]' ⊡ fuu
        Gux = B'*P[k+1]*A     + p[k+1]' ⊡ fux
        Gxu = copy(transpose(Gux))
    else
        Gxx = Q + A'*P[k+1]*A #+ p[k+1]' ⊡ fxx
        Guu = R + B'*P[k+1]*B #+ p[k+1]' ⊡ fuu
        Gux = B'*P[k+1]*A #+ p[k+1]' ⊡ fux
        Gxu = copy(transpose(Gux))
    end

    ρ = 0
    d[k] = (Guu+ρ*I)\gu
    K[k] = (Guu+ρ*I)\Gux

    p[k] = gx - K[k]'*gu + K[k]'*Guu*d[k] - Gxu*d[k]
    P[k] = Gxx + K[k]'*Guu*K[k] - Gxu*K[k] - K[k]'*Gux

    # ΔJ += gu'*d[k]
end

return P

end
