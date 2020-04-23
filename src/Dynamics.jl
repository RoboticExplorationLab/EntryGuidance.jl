function vinh_dynamics(x,u,p::PlanetModel)
    #Taken from pages 2-11 to 2-12 of Hypersonic Flight Mechanics by Busemann, Vinh, and Culp
    #https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19760024112.pdf

    #unpack state vector
    r = x[1]
    θ = x[2] #longitude
    ϕ = x[3] #latitude
    v = x[4] #velocity magnitude
    γ = x[5] #flight path angle
    ψ = x[6] #heading

    #unpack control vector
    D = u[1] #Drag acceleration magnitude
    L = u[2] #Lift acceleration magnitude
    σ = u[3] #Bank angle (specifies lift direction)

    #get gravity (assumed spherical for now)
    g = norm(gravitational_acceleration(p, [r,0.0,0.0]))

    #Terms involving Ω (planet rotation) are often thrown out.
    Ω = p.Ω #Set Ω = 0.0 here if you want that behavior

    ṙ = v*sin(γ)
    θ̇ = v*cos(γ)*cos(ψ)/(r*cos(ϕ))
    ϕ̇ = v*cos(γ)*sin(ψ)/r

    v̇ = -D - g*sin(γ) + Ω*Ω*r*cos(ϕ)*(sin(γ)*cos(ϕ) - cos(γ)*sin(ϕ)*sin(ψ))
    γ̇ = (L*cos(σ) - g*cos(γ) + (v*v/r)*cos(γ) + 2*Ω*v*cos(ϕ)*cos(ψ) + Ω*Ω*r*cos(ϕ)*(cos(γ)*cos(ϕ)+sin(γ)*sin(ϕ)*sin(ψ)))/v
    ψ̇ = (L*sin(σ)/cos(γ) - (v*v/r)*cos(γ)*cos(ψ)*tan(ϕ) + 2*Ω*v*(tan(γ)*cos(ϕ)*sin(ψ) - sin(ϕ)) - (Ω*Ω*r/cos(γ))*sin(ϕ)*cos(ϕ)*cos(ψ))/v

    ẋ = [ṙ, θ̇, ϕ̇, v̇, γ̇, ψ̇]
end

function cartesian_dynamics(x,u)

end

function quasilinear_dynamics(x,u)

end

# function full_dynamics(x,u)
#
# end
