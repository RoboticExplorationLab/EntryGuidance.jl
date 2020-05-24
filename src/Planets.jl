#Units are L=km M=kg T=hours

struct Planet{T}
    R::T #radius
    Ω::T #rotation rate (rad/T)
    gravity::AbstractGravityField
    atmosphere::AbstractAtmosphere
end

function Earth()
    p = Planet{Float64}(6378.1, 7.292115e-5*3600.0, EarthGravity(), EarthExponentialAtmosphere())
end

function EarthJ2()
    p = Planet{Float64}(6378.1, 7.292115e-5*3600.0, EarthGravityJ2(), EarthExponentialAtmosphere())
end

function Mars()
    p = Planet{Float64}(3396.2, 7.08824e-5*3600.0, MarsGravity(), MarsBiExponentialAtmosphere())
end

function MarsJ2()
    p = Planet{Float64}(3396.2, 7.08824e-5*3600.0, MarsaGravityJ2(), MarsBiExponentialAtmosphere())
end

function atmospheric_density(r::AbstractVector{T}, p::Planet{S}) where {T,S}
    atmospheric_density(r, p.atmosphere)
end

function atmospheric_density(r::T, p::Planet{S}) where {T,S}
    atmospheric_density(r, p.atmosphere)
end

function gravitational_acceleration(r::AbstractVector{T}, p::Planet{S}) where {T,S}
    gravitational_acceleration(r, p.gravity)
end
