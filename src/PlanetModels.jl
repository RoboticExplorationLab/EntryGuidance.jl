include("GravityFields.jl")
include("Atmospheres.jl")

struct PlanetModel{T}
    gravity::AbstractGravityField{T}
    atmosphere::AbstractAtmosphere{T}
end

function atmospheric_density(p::PlanetModel{T}, r::AbstractArray{T}) where {T}
    atmospheric_density(p.atmosphere, r)
end

function atmospheric_density(p::PlanetModel{T}, r::T) where {T}
    atmospheric_density(p.atmosphere, r)
end

function gravitation_acceleration(p::PlanetModel{T}, r::AbstractArray{T}) where {T}
    gravitation_acceleration(p.gravity, r)
end
