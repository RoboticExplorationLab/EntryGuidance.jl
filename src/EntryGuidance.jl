module EntryGuidance

include("GravityFields.jl")
export AbstractGravityField, SphericalGravityField, J2GravityField
export EarthGravity, EarthGravityJ2, MarsGravity, MarsGravityJ2
export gravitation_acceleration

include("Atmospheres.jl")
export AbstractAtmosphere, ExponentialAtmosphere
export EarthExponentialAtmosphere, MarsExponentialAtmosphere
export atmospheric_density

include("PlanetModels.jl")
export PlanetModel, EarthPlanetModel, EarthPlanetModelJ2, MarsPlanetModel, MarsPlanetModelJ2

include("Vehicles.jl")
export VehicleModel
export SimpleSphereConeVehicle

include("CoordinateTransformations.jl")
export planet_fixed_to_inertial, inertial_to_planet_fixed, cartesian_to_vinh, vinh_to_cartesian

include("Dynamics.jl")
export vinh_dynamics!, cartesian_dynamics!, bilinear_dynamics!

include("Controllers.jl")
export cartesian_no_lift, vinh_no_lift

end # module
