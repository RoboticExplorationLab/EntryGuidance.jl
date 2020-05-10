#Units are L=km M=kg T=hours

struct VehicleModel{T}
    m::T #mass
    A::T #area
end

function SimpleSphereConeVehicle()
    s = VehicleModel{Float64}(600.0, π*1.7e-3*1.7e-3)
end

function drag_coefficient(α::Number, s::VehicleModel)
    #This is valid for a 70 deg. sphere-cone above about Mach 8
    #From Gallais, "Atmospheric Re-Entry Vehicle Mechanics," p.82
    Cd = 1.65 - 2.77*α*α
end

function lift_coefficient(α::Number, s::VehicleModel)
    #This is valid for a 70 deg. sphere-cone over a wide Mach range
    #From Gallais, "Atmospheric Re-Entry Vehicle Mechanics," p.83
    Cl = 1.42*α
end

function drag_coefficient(s::VehicleModel)
    #This is super dumb for now
    Cd = 1.6
end

function lift_coefficient(s::VehicleModel)
    #This is super dumb for now
    Cl = 0.2
end
