module EntryGuidance

include("GravityFields.jl")
export AbstractGravityField, SphericalGravityField, J2GravityField
export EarthGravity, EarthGravityJ2, MarsGravity, MarsGravityJ2
export gravitational_acceleration

include("Atmospheres.jl")
export AbstractAtmosphere, ExponentialAtmosphere, BiExponentialAtmosphere
export EarthExponentialAtmosphere, MarsExponentialAtmosphere, MarsBiExponentialAtmosphere
export atmospheric_density

include("Planets.jl")
export Planet, Earth, EarthJ2, Mars, MarsJ2
export gravitational_acceleration, atmospheric_density

include("Vehicles.jl")
export AbstractVehicle, HypersonicVehicle
export SimpleMSLVehicle
export drag_coefficient, lift_coefficient

include("Models.jl")
export AbstractModel, VinhModel, CartesianModel
export CartesianMSLModel, CartesianMSLModelJ2, VinhMSLModel, VinhMSLModelJ2
export drag_coefficient, lift_coefficient, atmospheric_density, gravitational_acceleration

include("CoordinateTransformations.jl")
export planet_fixed_to_inertial, inertial_to_planet_fixed, cartesian_to_vinh, vinh_to_cartesian

include("Dynamics.jl")
export dynamics!, angles_input

include("Controllers.jl")


end # module
