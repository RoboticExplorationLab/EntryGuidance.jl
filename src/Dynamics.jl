using LinearAlgebra

function dynamics!(ẋ::AbstractVector{T},x::AbstractVector{T},u::AbstractVector{T},m::VinhModel{T}) where {T}
    #Taken from pages 2-11 to 2-12 of Hypersonic Flight Mechanics by Busemann, Vinh, and Culp
    #https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19760024112.pdf

    #unpack state vector
    r = x[1]
    θ = x[2] #longitude
    ϕ = x[3] #latitude
    v = x[4] #velocity magnitude
    γ = x[5] #flight path angle
    ψ = x[6] #heading

    #unpack control vector
    D = u[1] #Drag acceleration magnitude
    L = u[2] #Lift acceleration magnitude
    σ = u[3] #Bank angle (specifies lift direction)

    #get gravity (assumed spherical for now)
    g = norm(gravitational_acceleration([r,0.0,0.0], m))

    #Terms involving Ω (planet rotation) are often thrown out in the literature.
    Ω = m.planet.Ω #Set Ω = 0.0 here if you want that behavior

    ṙ = v*sin(γ)
    θ̇ = v*cos(γ)*cos(ψ)/(r*cos(ϕ))
    ϕ̇ = v*cos(γ)*sin(ψ)/r

    v̇ = -D - g*sin(γ) + Ω*Ω*r*cos(ϕ)*(sin(γ)*cos(ϕ) - cos(γ)*sin(ϕ)*sin(ψ))
    γ̇ = (L*cos(σ) - g*cos(γ) + (v*v/r)*cos(γ) + 2*Ω*v*cos(ϕ)*cos(ψ) + Ω*Ω*r*cos(ϕ)*(cos(γ)*cos(ϕ)+sin(γ)*sin(ϕ)*sin(ψ)))/v
    ψ̇ = (L*sin(σ)/cos(γ) - (v*v/r)*cos(γ)*cos(ψ)*tan(ϕ) + 2*Ω*v*(tan(γ)*cos(ϕ)*sin(ψ) - sin(ϕ)) - (Ω*Ω*r/cos(γ))*sin(ϕ)*cos(ϕ)*cos(ψ))/v

    ẋ .= [ṙ, θ̇, ϕ̇, v̇, γ̇, ψ̇]
end

function dynamics!(ẋ::AbstractVector{T},x::AbstractVector{T},u::AbstractVector{T},m::CartesianModel{T}) where {T}

    #unpack state vector
    r = x[1:3]
    v = x[4:6]

    #unpack control vector
    D = u[1] #Drag acceleration magnitude
    L = u[2] #Lift acceleration magnitude
    σ = u[3] #Bank angle (specifies lift direction)

    #get gravity
    g = gravitational_acceleration(r, m)

    #Terms involving Ω (planet rotation) are often thrown out in the literature.
    Ω = m.planet.Ω #Set Ω = 0.0 here if you want that behavior
    Ω̂ = hat([0, 0, Ω])

    #Aerodynamic acceleration
    #Bank angle is computed relative to "vertical" (r, v) plane.
    #σ=0 corresponds to lift poining in the "near radial" direction
    #Positive σ rotates lift vector around v in right-hand sense towards r×v
    e1 = cross(r,v)
    e1 = e1./norm(e1)
    e2 = cross(v,e1)
    e2 = e2./norm(e2)
    a = -(D/norm(v)).*v + L*sin(σ)*e1 + L*cos(σ)*e2

    v̇ = a + g - 2*Ω̂*v - Ω̂*Ω̂*r

    ẋ .= [v; v̇]
end

function angles_input(u::AbstractVector{T},x::AbstractVector{T},model::VinhModel{T}) where {T}

    #unpack control
    α = u[1] #angle of attack
    σ = u[2] #bank angle

    #unpack state
    r = x[1]
    v = x[4]

    #Assume spherically symmetric atmosphere
    ρ = atmospheric_density([r,0,0], model)

    #Calculate drag acceleration
    Cd = drag_coefficient(α, model)
    A = model.vehicle.A
    m = model.vehicle.m
    D = 0.5*Cd*ρ*A*v*v/m

    #Calculate lift acceleration
    Cl = lift_coefficient(α, model)
    L = 0.5*Cl*ρ*A*v*v/m

    a = [D,L,σ]
end

function angles_input(u::AbstractVector{T},x::AbstractVector{T},model::CartesianModel{T}) where {T}

    #unpack control
    α = u[1] #angle of attack
    σ = u[2] #bank angle

    #unpack state
    r = x[1:3]
    v = x[4:6]
    V = norm(v)

    #atmospheric density
    ρ = atmospheric_density(r, model)

    #Calculate drag acceleration
    Cd = drag_coefficient(α, model)
    A = model.vehicle.A
    m = model.vehicle.m
    D = 0.5*Cd*ρ*A*V*V/m

    #Calculate lift acceleration
    Cl = lift_coefficient(α, model)
    L = 0.5*Cl*ρ*A*V*V/m

    a = [D,L,σ]
end

# function linear_dynamics!(ẋ::AbstractVector{T},x::AbstractVector{T},u::AbstractVector{T},y::AbstractVector{T},p::PlanetModel{T}) where {T}
#
#     #unpack state vector
#     r = x[1:3]
#     v = x[4:6]
#
#     #unpack parameter vector
#     kg = y[1] # = g0/r0 (equivalent spring force for gravity - assumes variations in g and r are small)
#
#     #Terms involving Ω (planet rotation) are often thrown out in the literature.
#     Ω = p.Ω #Set Ω = 0.0 here if you want that behavior
#     Ω̂ = hat([0, 0, Ω])
#
#     A = [zeros(3,3) I; -kg*I-Ω̂*Ω̂ -2*Ω̂]
#     B = [zeros(3,3); I]
#
#     ẋ .= A*x + B*u
# end
#
# function quasilinear_dynamics!(ẋ::AbstractVector{T},x::AbstractVector{T},u::AbstractVector{T},p::PlanetModel{T}) where {T}
#
#     #unpack state vector
#     r = x[1:3]
#     v = x[4:6]
#
#     #calculate gravity term
#     kg = norm(gravitational_acceleration(r,p))/norm(r) #equivalent spring force for gravity - assumes variations in g and r are small
#
#     #Terms involving Ω (planet rotation) are often thrown out in the literature.
#     Ω = p.Ω #Set Ω = 0.0 here if you want that behavior
#     Ω̂ = hat([0, 0, Ω])
#
#     A = [zeros(3,3) I; -kg*I-Ω̂*Ω̂ -2*Ω̂]
#     B = [zeros(3,3); I]
#
#     ẋ .= A*x + B*u
# end
