#Units are L=km M=kg T=hours

abstract type AbstractGravityField{T} end

struct SphericalGravityField{T} <: AbstractGravityField{T}
    μ::T #standard gravitational parameter
end

struct J2GravityField{T} <: AbstractGravityField{T}
    μ::T #standard gravitational parameter (L^3/T^2)
    R::T #planet radius (L)
    J2::T #J2 coefficient (dimensionless)
end

function EarthGravity()
    SphericalGravityField{Float64}(398600.4418*(3600.0^2))
end

function MarsGravity()
    SphericalGravityField{Float64}(42828.37*(3600.0^2))
end

function EarthGravityJ2()
    J2GravityField{FLoat64}(398600.4415*(3600.0^2), 6378.1363, 0.1082635854e-2)
end

function MarsGravityJ2()
    J2GravityField{Float64}(42828.37*(3600.0^2), 3396.203986, 0.196045e-2)
end

function gravitational_acceleration(r::AbstractVector{T}, g::SphericalGravityField{T}) where {T}
    R = norm(r)
    a = (-g.μ/(R*R*R)).*r
end

function gravitational_acceleration(r::AbstractVector{T}, g::J2GravityField{T}) where {T}
    R = norm(r)
    R3 = R^3
    R7 = R^7
    x = r[1]
    y = r[2]
    z = r[3]
    x2y2 = x*x + y*y
    z2 = z*z
    J2 = g.J2*(g.μ*g.R*g.R) #convert from dimensionless J2 to L^5/T^2 units
    a = Diagonal([-g.μ/R3 + J2*(6.0*z2 - 1.5*x2y2)/R7, -g.μ/R3 + J2*(6.0*z2 - 1.5*x2y2)/R7, -g.μ/R3 + J2*(3.0*z2 - 4.5*x2y2)/R7])*r
end
