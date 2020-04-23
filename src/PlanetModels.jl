
function EarthPlanetModel()
    p = PlanetModel{Float64}(EarthGravity(), EarthExponentialAtmosphere(), 7.292115e-5)
end

function EarthPlanetModelJ2()
    p = PlanetModel{Float64}(EarthGravityJ2(), EarthExponentialAtmosphere(), 7.292115e-5)
end

function MarsPlanetModel()
    p = PlanetModel{Float64}(MarsGravity(), MarsExponentialAtmosphere(), 7.08824e-5)
end

function MarsPlanetModelJ2()
    p = PlanetModel{Float64}(MarsaGravityJ2(), MarsExponentialAtmosphere(), 7.08824e-5)
end

struct PlanetModel{T}
    gravity::AbstractGravityField{T}
    atmosphere::AbstractAtmosphere{T}
    Ω::T #rotation rate (rad/sec)
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
