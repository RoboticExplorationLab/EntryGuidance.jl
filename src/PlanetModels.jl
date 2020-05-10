#Units are L=km M=kg T=hours

function EarthPlanetModel()
    p = PlanetModel{Float64}(6378.1, 7.292115e-5*3600.0, EarthGravity(), EarthExponentialAtmosphere())
end

function EarthPlanetModelJ2()
    p = PlanetModel{Float64}(6378.1, 7.292115e-5*3600.0, EarthGravityJ2(), EarthExponentialAtmosphere())
end

function MarsPlanetModel()
    p = PlanetModel{Float64}(3396.2, 7.08824e-5*3600.0, MarsGravity(), MarsBiExponentialAtmosphere())
end

function MarsPlanetModelJ2()
    p = PlanetModel{Float64}(3396.2, 7.08824e-5*3600.0, MarsaGravityJ2(), MarsBiExponentialAtmosphere())
end

struct PlanetModel{T}
    R::T #radius
    Ω::T #rotation rate (rad/T)
    gravity::AbstractGravityField{T}
    atmosphere::AbstractAtmosphere{T}
end

function atmospheric_density(r::AbstractVector{T}, p::PlanetModel{T}) where {T}
    atmospheric_density(r, p.atmosphere)
end

function atmospheric_density(r::T, p::PlanetModel{T}) where {T}
    atmospheric_density(r, p.atmosphere)
end

function gravitational_acceleration(r::AbstractVector{T}, p::PlanetModel{T}) where {T}
    gravitational_acceleration(r, p.gravity)
end
