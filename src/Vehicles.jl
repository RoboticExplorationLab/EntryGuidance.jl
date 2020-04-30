
function SimpleSphereConeVehicle()
    #kg, km^2
    s = VehicleModel{Float64}(600.0, π*1.7e-3*1.7e-3)
end

struct VehicleModel{T}
    m::T #mass
    A::T #area
end

function drag_coefficient(s::VehicleModel)
    #This is super dumb for now
    Cd = 2.0
end

function lift_coefficient(s::VehicleModel)
    #This is super dumb for now
    Cl = 0.2
end
