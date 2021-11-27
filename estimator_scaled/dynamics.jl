struct EntryVehicle{T}
    evmodel::EG.CartesianModel{T}
    dscale::Float64
    tscale::Float64
    uscale::Float64
end
function my_atmo(model,x,H)
    r = norm(x[1:3])
    r0 =3396.2
    ρ0 = 5.25*1e7
    # H =11.1
    # ρ0 *= 1e7
    ρ = ρ0*exp(-(r-r0)/H)
    return ρ
end
function evdynamics(model::EntryVehicle, x, u)

    # unpack scaling
    dscale, tscale, uscale = model.dscale, model.tscale, model.uscale

    r = x[SA[1,2,3]]*dscale
    v = x[SA[4,5,6]]*(dscale/tscale)
    σ = x[7]
    kρ = x[8]
    # ρ0 = x[9]
    α = deg2rad(15)
    # σ̇ = u
    #atmospheric density
    # ρ2 = atmospheric_density(r, model.evmodel)
    # ρ = atmospheric_density(r,model.exp_at)
    ρ = my_atmo(model, r, 7.295)*kρ

    # @show (ρ - ρ2)/ρ2
    # @show ρ2
    # @show ρ

    #Calculate drag acceleration
    Cd = drag_coefficient(α, model.evmodel)
    A = model.evmodel.vehicle.A
    m = model.evmodel.vehicle.m
    D = 0.5*Cd*ρ*A*dot(v,v)/m

    #Calculate lift acceleration
    Cl = lift_coefficient(α, model.evmodel)
    L = 0.5*Cl*ρ*A*dot(v,v)/m

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
    # L_a = e1*u[1] + e2*u[2]
    L_a = L*sin(σ)*e1 + L*cos(σ)*e2
                      # this is rotating planet effects
    v̇ = D_a + L_a + g - 2*Ω̂*v - Ω̂*Ω̂*r

    # # epsilon dot stuff
    # mu =  model.evmodel.planet.gravity.μ
    # vi = v + Ω̂*r
    # dvi = v̇ + Ω̂*v
    # epsilon_dot = ( dot(dvi,vi) + mu*dot(v,r)/norm(r)^(3) )
    # epsilon_dot /= 1e8 # scaling

    # rescale units
    v = v/(dscale/tscale)
    v̇ = v̇/(dscale/tscale^2)

    return SA[v[1],v[2],v[3],v̇[1],v̇[2],v̇[3],u[1]*1000,0]
end

function rk4(model,x_n,u,dt)
    k1 = dt*evdynamics(model,x_n,u)
    k2 = dt*evdynamics(model,x_n+k1/2,u)
    k3 = dt*evdynamics(model,x_n+k2/2,u)
    k4 = dt*evdynamics(model,x_n+k3,u)
    return (x_n + (1/6)*(k1 + 2*k2 + 2*k3 + k4))
end
function getAB(model,X,U,dt)
    # N = length(X)
    # n = length(X[1])
    # m = length(U[1])
    # A = [@SArray zeros(n,n) for i = 1:(N-1)]
    # B = [@SArray zeros(n,m) for i = 1:(N-1)]
    # for k = 1:(N-1)
        A= ForwardDiff.jacobian(_x -> rk4(model,_x,U,dt),X)
        B= ForwardDiff.jacobian(_u -> rk4(model,X,_u,dt),U)
    # end
    return A,B
end
