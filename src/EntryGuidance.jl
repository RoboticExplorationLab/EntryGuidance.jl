module EntryGuidance

include("GravityFields.jl")
include("Atmospheres.jl")
include("PlanetModels.jl")
include("CoordinateTransformations.jl")
include("Dynamics.jl")
export AbstractGravityField, SphericalGravityField, J2GravityField
export EarthGravity, EarthGravityJ2, MarsGravity, MarsGravityJ2
export AbstractAtmosphere, ExponentialAtmosphere
export EarthExponentialAtmosphere, MarsExponentialAtmosphere
export PlanetModel, EarthPlanetModel, EarthPlanetModelJ2, MarsPlanetModel, MarsPlanetModelJ2
export gravitation_acceleration, atmospheric_density, vinh_dynamics, cartesian_dynamics, quasilinear_dynamics

end # module
