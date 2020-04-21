module EntryGuidance

include("PlanetModels.jl")
export AbstractGravityField, SphericalGravityField, J2GravityField
export EarthGravity, EarthGravityJ2, MarsGravity, MarsGravityJ2
export AbstractAtmosphere, ExponentialAtmosphere
export PlanetModel

end # module
