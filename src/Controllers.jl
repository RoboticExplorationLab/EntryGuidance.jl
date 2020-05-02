using LinearAlgebra

function vinh_fixed_lift(x::AbstractVector{T},s::VehicleModel{T},p::PlanetModel{T}) where {T}
    #Open-loop with L=0

    #unpack state
    r = x[1]
    v = x[4]

    #Assume spherically symmetric atmosphere
    ρ = atmospheric_density([r,0,0], p)

    #Calculate drag acceleration
    Cd = drag_coefficient(s)
    A = s.A
    m = s.m
    D = 0.5*Cd*ρ*A*v*v/m

    #Calculate lift acceleration
    Cl = lift_coefficient(s)
    L = 0.5*Cl*ρ*A*v*v/m
    σ = π/4

    u = [D,L,σ]
end

function cartesian_fixed_lift(x::AbstractVector{T},s::VehicleModel{T},p::PlanetModel{T}) where {T}
    #Open-loop with L=0

    #unpack state
    r = x[1:3]
    v = x[4:6]
    V = norm(v)

    #atmospheric density
    ρ = atmospheric_density(r, p)

    #Calculate drag acceleration
    Cd = drag_coefficient(s)
    A = s.A
    m = s.m
    D = 0.5*Cd*ρ*A*V*V/m

    #Calculate lift acceleration
    Cl = lift_coefficient(s)
    L = 0.5*Cl*ρ*A*V*V/m
    σ = π/4

    u = [D,L,σ]
end
