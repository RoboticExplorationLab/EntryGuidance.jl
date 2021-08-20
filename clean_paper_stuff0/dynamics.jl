function getmaxL(model,x)
    # gets the maximum available lift magnitude given position and velocity

    # unpack state
    r = x[1:3]
    v = x[4:6]

    # get density
    ρ = atmospheric_density(r, model.evmodel)
    A = model.evmodel.vehicle.A
    m = model.evmodel.vehicle.m

    # calculate lift
    Cl = lift_coefficient(deg2rad(20), model.evmodel)
    L = 0.5*Cl*ρ*A*dot(v,v)/m

    # return scaled version
    return L/model.uscale
end
struct EntryVehicle{T}
    evmodel::EG.CartesianModel{T}
    uscale::Float64
end
function evdynamics(model::EntryVehicle, x, u)

    #unpack state
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
    Ω = model.evmodel.planet.Ω #Set Ω = 0.0 here if you want that behavior
    Ω̂ = hat([0, 0, Ω])

    #Aerodynamic acceleration
    e1 = cross(r,v)
    e1 = e1/norm(e1)
    e2 = cross(v,e1)
    e2 = e2/norm(e2)
    D_a = -(D/norm(v))*v #+ L*sin(σ)*e1 + L*cos(σ)*e2
    L_mag = getmaxL(model,x)
    L_a = L_mag*(e1*u[1] + e2*u[2])
                      # this is rotating planet effects
    v̇ = D_a + L_a + g - 2*Ω̂*v - Ω̂*Ω̂*r

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
