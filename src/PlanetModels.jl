struct PlanetModel{T}
    gravity::AbstractGravityField{T}
    atmosphere::AbstractAtmosphere{T}
end

abstract type AbstractGravityField{T} end

struct SphericalGravityField{T} <: AbstractGravityField{T}
    μ::T #standard gravitational parameter
end

struct J2GravityField{T} <: AbstractGravityField{T}
    μ::T #standard gravitational parameter
    J2::T #J2 coefficient
end

function gravitational_acceleration(g::SphericalGravityField{T}, r::AbstractArray{T}) where {T}
    R = norm(r)
    a = (-g.μ/(R*R*R)).*r
end

function gravitational_acceleration(g::J2GravityField{T}, r::AbstractArray{T}) where {T}
    R = norm(r)
    R3 = R^3
    R7 = R^7
    x = r[1]
    y = r[2]
    z = r[3]
    x2y2 = x*x + y*y
    z2 = z*z
    a = Diagonal([-g.μ/R3 + g.J2*(6.0*z2 - 1.5*x2y2)/R7, -g.μ/R3 + g.J2*(6.0*z2 - 1.5*x2y2)/R7, -g.μ/R3 + g.J2*(3.0*z2 - 4.5*x2y2)/R7])*r
end

function gravitation_acceleration(p::PlanetModel{T}, r::AbstractArray{T}) where {T}
    gravitation_acceleration(p.gravity, r)
end

abstract type AbstractAtmosphere{T} end

struct ExponentialAtmosphere{T} <: AbstractAtmosphere{T}
    r0::T #planet radius
    H::T #scale height
    ρ0::T #surface density
end

function atmospheric_density(a::ExponentialAtmosphere{T}, r::AbstractArray{T}) where {T}
    R = norm(r)
    ρ = a.ρ0*exp(-(R-a.r0)/a.H)
end

function atmospheric_density(a::ExponentialAtmosphere{T}, r::T) where {T}
    ρ = a.ρ0*exp(-(r-a.r0)/a.H)
end

function atmospheric_density(p::PlanetModel{T}, r::AbstractArray{T}) where {T}
    atmospheric_density(p.atmosphere, r)
end

function atmospheric_density(p::PlanetModel{T}, r::T) where {T}
    atmospheric_density(p.atmosphere, r)
end
