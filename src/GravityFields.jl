abstract type AbstractGravityField{T} end

# function EarthGravity()
#     #Units are Kg, Km, s
#     SphericalGravityField{Float64}(398600.4418)
# end

function MarsGravity()
    #Units are Kg, Km, s
    SphericalGravityField{Float64}(42828.37)
end

function EarthGravityJ2()
    #Units are Kg, Km, s
    J2GravityField{FLoat64}(398600.4415, 6378.1363, 0.1082635854e-2)
end

function MarsGravityJ2()
    #Units are Kg, Km, s
    J2GravityField{Float64}(42828.37, 3396.0, 0.196045e-2)
end

struct SphericalGravityField{T} <: AbstractGravityField{T}
    μ::T #standard gravitational parameter
end

struct J2GravityField{T} <: AbstractGravityField{T}
    μ::T #standard gravitational parameter ()
    R::T #planet radius (L)
    J2::T #dimensionless J2 coefficient
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
    J2 = g.J2*(g.μ*g.R*g.R) #convert from dimensionless J2 to L^5/T^2 units
    a = Diagonal([-g.μ/R3 + J2*(6.0*z2 - 1.5*x2y2)/R7, -g.μ/R3 + J2*(6.0*z2 - 1.5*x2y2)/R7, -g.μ/R3 + J2*(3.0*z2 - 4.5*x2y2)/R7])*r
end
