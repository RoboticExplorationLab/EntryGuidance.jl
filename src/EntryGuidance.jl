module EntryGuidance

include("GravityFields.jl")
export AbstractGravityField, SphericalGravityField, J2GravityField
export EarthGravity, EarthGravityJ2, MarsGravity, MarsGravityJ2
export gravitational_acceleration

include("Atmospheres.jl")
export AbstractAtmosphere, ExponentialAtmosphere, BiExponentialAtmosphere
export EarthExponentialAtmosphere, MarsExponentialAtmosphere, MarsBiExponentialAtmosphere
export atmospheric_density

include("PlanetModels.jl")
export PlanetModel, EarthPlanetModel, EarthPlanetModelJ2, MarsPlanetModel, MarsPlanetModelJ2
export gravitational_acceleration, atmospheric_density

include("Vehicles.jl")
export VehicleModel
export SimpleSphereConeVehicle
export drag_coefficient, lift_coefficient

include("CoordinateTransformations.jl")
export planet_fixed_to_inertial, inertial_to_planet_fixed, cartesian_to_vinh, vinh_to_cartesian

include("Dynamics.jl")
export vinh_dynamics!, cartesian_dynamics!, linear_dynamics!, quasilinear_dynamics!

include("Controllers.jl")
export cartesian_test_controller, vinh_test_controller, linear_test_controller

end # module
