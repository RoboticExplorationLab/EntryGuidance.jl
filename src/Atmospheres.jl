abstract type AbstractAtmosphere{T} end

function EarthExponentialAtmosphere()
    #Units are Kg, Km, s
    #https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html
    ExponentialAtmosphere{Float64}(6378.1, 8.5, 1.217e9)
end

function MarsExponentialAtmosphere()
    #Units are Km, Kg, s
    #https://nssdc.gsfc.nasa.gov/planetary/factsheet/marsfact.html
    ExponentialAtmosphere{Float64}(3396.2, 11.1, 2.0e7)
end

struct ExponentialAtmosphere{T} <: AbstractAtmosphere{T}
    r0::T #reference radius (L)
    H::T #scale height (L)
    ρ0::T #surface density (M/L^3)
end

function atmospheric_density(a::ExponentialAtmosphere{T}, r::AbstractArray{T}) where {T}
    R = norm(r)
    ρ = a.ρ0*exp(-(R-a.r0)/a.H)
end

function atmospheric_density(a::ExponentialAtmosphere{T}, r::T) where {T}
    ρ = a.ρ0*exp(-(r-a.r0)/a.H)
end
