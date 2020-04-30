using LinearAlgebra

function vinh_no_lift(x::AbstractVector{T},s::VehicleModel{T},p::PlanetModel{T}) where {T}
    #Open-loop with L=0

    #grab velocity magnitude from state
    v = x[4]

    #Assume spherically symmetric atmosphere
    ρ = atmospheric_density([r,0,0], p)

    #Calculate drag acceleration
    Cd = drag_coefficient(s)
    A = s.A
    m = s.m
    D = 0.5*Cd*ρ*A*v*v/m

    u = [D,0,0]
end

function cartesian_no_lift(x::AbstractVector{T},s::VehicleModel{T},p::PlanetModel{T}) where {T}
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

    u = [D,0,0]
end
