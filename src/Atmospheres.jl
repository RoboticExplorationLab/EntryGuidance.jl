#Units are L=km M=kg T=hours

abstract type AbstractAtmosphere end

struct ExponentialAtmosphere{T} <: AbstractAtmosphere
    r0::T #reference radius (L)
    H::T #scale height (L)
    ρ0::T #surface density (M/L^3)
end

struct BiExponentialAtmosphere{T} <: AbstractAtmosphere
    r0::T #reference radius (L)
    transition_altitude::T #transition altitude (L)
    H1::T #low altitude scale height (L)
    ρ1::T #low altitude reference density (M/L^3)
    H2::T #high altitude scale height (L)
    ρ2::T #high altitude reference density (M/L^3)
end

function EarthExponentialAtmosphere()
    #https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html
    ExponentialAtmosphere{Float64}(6378.1, 8.5, 1.217e9)
end

function MarsExponentialAtmosphere()
    #https://nssdc.gsfc.nasa.gov/planetary/factsheet/marsfact.html
    ExponentialAtmosphere{Float64}(3396.2, 11.1, 2.0e7)
end

function MarsBiExponentialAtmosphere()
    #From Gallais, "Atmospheric Re-Entry Vehicle Mechanics," p.29
    BiExponentialAtmosphere{Float64}(3396.2, 25.0, 11.049, 1.59e7, 7.295, 5.25e7)
end

function atmospheric_density(r::AbstractVector{T}, a::ExponentialAtmosphere{S}) where {T,S}
    R = norm(r)
    atmospheric_density(R,a)
end

function atmospheric_density(r::T, a::ExponentialAtmosphere{S}) where {T,S}
    ρ = a.ρ0*exp(-(r-a.r0)/a.H)
end

function atmospheric_density(r::AbstractVector{T}, a::BiExponentialAtmosphere{S}) where {T,S}
    R = norm(r)
    atmospheric_density(R,a)
end

function atmospheric_density(r::T, a::BiExponentialAtmosphere{S}) where {T,S}
    #Smoothly interpolate between two exponential fits
    alt = r-a.r0
    ρ_low = a.ρ1*exp(-alt/a.H1)
    ρ_high = a.ρ2*exp(-alt/a.H2)
    β = 0.5*(1.0+tanh((alt-a.transition_altitude)))
    ρ = (1-β)*ρ_low + β*ρ_high
end
