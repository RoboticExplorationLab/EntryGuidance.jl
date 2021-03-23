
function ezdynamics(x, u)

    #unpack state
    r = x[1:3]
    v = x[4:6]

    #atmospheric density
    ρ = 1.2*10 # kg/m³

    #Calculate drag acceleration
    Cd = 2.2
    A = 1.0  # m²
    m = 1000 # kg
    D = 0.5*Cd*ρ*A*dot(v,v)/m

    # this is maximum allowable lift
    # Cl = 1.42*α (in radians)
    # L = 0.5*Cl*ρ*A*V*V/m

    #get gravity
    g = [0;0;-9.8]

    #Aerodynamic acceleration
    e1 = cross(r,v)
    e1 .= e1/norm(e1)
    e2 = cross(v,e1)
    e2 .= e2/norm(e2)

    # lift vector (we get to control this within allowable bounds)
    L_a = e1*u[1] + e2*u[2]
    D_a = -(D/norm(v))*v

    v̇ = L_a + D_a + g #- 2*Ω̂*v - Ω̂*Ω̂*r

    return [v; v̇]
end
function rk4(x_n,u,dt)
    k1 = dt*ezdynamics(x_n,u)
    k2 = dt*ezdynamics(x_n+k1/2,u)
    k3 = dt*ezdynamics(x_n+k2/2,u)
    k4 = dt*ezdynamics(x_n+k3,u)
    return (x_n + (1/6)*(k1+ 2*k2 + 2*k3 + k4))
end

function getAB(X,U,dt)
    N = length(X)
    n = length(X[1])
    m = length(U[1])
    A = [zeros(n,n) for i = 1:(N-1)]
    B = [zeros(n,m) for i = 1:(N-1)]
    for k = 1:(N-1)
        A[k] = ForwardDiff.jacobian(_x -> rk4(_x,U[k],dt),X[k])
        B[k] = ForwardDiff.jacobian(_u -> rk4(X[k],_u,dt),U[k])
    end
    return A,B
end
