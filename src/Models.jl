abstract type AbstractModel{T} end

struct VinhModel{T} <: AbstractModel{T}
    planet::Planet{T}
    vehicle::AbstractVehicle{T}
end

struct CartesianModel{T} <: AbstractModel{T}
    planet::Planet{T}
    vehicle::AbstractVehicle{T}
end

function CartesianMSLModel()
    m = CartesianModel(Mars(), SimpleMSLVehicle())
end

function CartesianMSLModelJ2()
    m = CartesianModel(MarsJ2(), SimpleMSLVehicle())
end

function VinhMSLModel()
    m = VinhModel(Mars(), SimpleMSLVehicle())
end

function VinhMSLModelJ2()
    m = VinhModel(MarsJ2(), SimpleMSLVehicle())
end

function drag_coefficient(α::T, m::AbstractModel{T}) where {T}
    #Quadratic in angle of attack (valid up to ~10 deg.)
    Cd = drag_coefficient(α, m.vehicle)
end

function lift_coefficient(α::T, m::AbstractModel{T}) where {T}
    #Linear in angle of attack (valid up to ~10 deg.)
    Cl = lift_coefficient(α, m.vehicle)
end

function atmospheric_density(r::AbstractVector{T}, m::AbstractModel{T}) where {T}
    atmospheric_density(r, m.planet)
end

function atmospheric_density(r::T, m::AbstractModel{T}) where {T}
    atmospheric_density(r, m.planet)
end

function gravitational_acceleration(r::AbstractVector{T}, m::AbstractModel{T}) where {T}
    gravitational_acceleration(r, m.planet)
end
