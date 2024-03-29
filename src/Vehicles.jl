#Units are L=km M=kg T=hours

abstract type AbstractVehicle end

struct HypersonicVehicle{T} <: AbstractVehicle
    m::T #mass (M)
    A::T #area (L^2)
    Cd0::T #drag at zero angle of attack
    Cd2::T #quadratic drag coefficient
    Cl1::T #linear lift coefficient
end

function SimpleMSLVehicle()
    #Mass and area are from MSL
    #Aero coefficients are valid for a hypersonic 70 deg. sphere-cone
    #From Gallais, "Atmospheric Re-Entry Vehicle Mechanics," p.82-83
    s = HypersonicVehicle{Float64}(2400.0, π*2.25e-3*2.25e-3, 1.65, -2.77, 1.42)
end

function drag_coefficient(α::T, s::HypersonicVehicle{S}) where {T,S}
    #Quadratic in angle of attack (valid up to ~10 deg.)
    Cd = s.Cd0 + s.Cd2*α*α
end

function lift_coefficient(α::T, s::HypersonicVehicle{S}) where {T,S}
    #Linear in angle of attack (valid up to ~10 deg.)
    Cl = s.Cl1*α
end
