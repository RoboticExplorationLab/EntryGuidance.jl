function vinh_dynamics!(ẋ::AbstractVector{T},x::AbstractVector{T},u::AbstractVector{T},p::PlanetModel{T}) where {T}
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
    g = norm(gravitational_acceleration([r,0.0,0.0], p))

    #Terms involving Ω (planet rotation) are often thrown out in the literature.
    Ω = p.Ω #Set Ω = 0.0 here if you want that behavior

    ṙ = v*sin(γ)
    θ̇ = v*cos(γ)*cos(ψ)/(r*cos(ϕ))
    ϕ̇ = v*cos(γ)*sin(ψ)/r

    v̇ = -D - g*sin(γ) + Ω*Ω*r*cos(ϕ)*(sin(γ)*cos(ϕ) - cos(γ)*sin(ϕ)*sin(ψ))
    γ̇ = (L*cos(σ) - g*cos(γ) + (v*v/r)*cos(γ) + 2*Ω*v*cos(ϕ)*cos(ψ) + Ω*Ω*r*cos(ϕ)*(cos(γ)*cos(ϕ)+sin(γ)*sin(ϕ)*sin(ψ)))/v
    ψ̇ = (L*sin(σ)/cos(γ) - (v*v/r)*cos(γ)*cos(ψ)*tan(ϕ) + 2*Ω*v*(tan(γ)*cos(ϕ)*sin(ψ) - sin(ϕ)) - (Ω*Ω*r/cos(γ))*sin(ϕ)*cos(ϕ)*cos(ψ))/v

    ẋ .= [ṙ, θ̇, ϕ̇, v̇, γ̇, ψ̇]
end

function cartesian_dynamics!(ẋ::AbstractVector{T},x::AbstractVector{T},u::AbstractVector{T},p::PlanetModel{T}) where {T}

    #unpack state vector
    r = x[1:3]
    v = x[4:6]

    #unpack control vector
    D = u[1] #Drag acceleration magnitude
    L = u[2] #Lift acceleration magnitude
    σ = u[3] #Bank angle (specifies lift direction)

    #get gravity
    g = gravitational_acceleration(r, p)

    #Terms involving Ω (planet rotation) are often thrown out in the literature.
    Ω = p.Ω #Set Ω = 0.0 here if you want that behavior
    Ω̂ = hat([0, 0, Ω])

    #Aerodynamic acceleration
    #Bank angle is computed relative to "vertical" (r, v) plane.
    #Positive σ rotates lift vector around v in right-hand sense towards r×v
    e1 = cross(r,v)
    e1 = e1./norm(e1)
    e2 = cross(e1,v)
    e2 = e2./norm(e2)
    a = -(D/norm(v)).*v + L*sin(σ)*e1 + L*cos(σ)*e2

    v̇ = a + g - 2*Ω̂*v - Ω̂*Ω̂*r

    ẋ .= [v; v̇]
end

function bilinear_dynamics!(ẋ::AbstractVector{T},x::AbstractVector{T},u::AbstractVector{T},y::AbstractVector{T},p::PlanetModel{T}) where {T}

    #unpack state vector
    r = x[1:3]
    v = x[4:6]

    #unpack parameter vector
    r0 = y[1]
    v0 = y[2]
    g0 = y[3]
    D0 = y[4]

    #Terms involving Ω (planet rotation) are often thrown out in the literature.
    Ω = p.Ω #Set Ω = 0.0 here if you want that behavior
    Ω̂ = hat([0, 0, Ω])

    #Bilinear dynamics model
    kg = g0/r0 #equivalent spring for gravity
    kd = D0/v0 #equivalent viscous damper for drag
    A = [zeros(3,3) I; -kg*I-Ω̂*Ω̂ -kd*I-2*Ω̂]
    B = [zeros(3,3); hat(v)] #bilinear term v×u

    ẋ .= A*x + B*u
end
