function hat(x)
    x̂ = [0 -x[3] x[2];
         x[3] 0 -x[1];
         -x[2] x[1] 0]
end

function planet_fixed_to_inertial(t::T,x::AbstractVector{T},p::Planet{T}) where {T}
    #Assume axes are aligned at t=0
    ω̂ = hat([0, 0, p.Ω])
    R = exp(t*ω̂)
    xn = [R*x[1:3]; R*(x[4:6] + ω̂*x[1:3])]
end

function planet_fixed_to_inertial(x::AbstractVector{T},p::Planet{T}) where {T}
    #Assume axes are aligned
    planet_fixed_to_inertial(0.0,x,p)
end

function inertial_to_planet_fixed(t::T,x::AbstractVector{T},p::Planet{T}) where {T}
    #Assume axes are aligned at t=0
    ω̂ = hat([0, 0, p.Ω])
    R = exp(-t*ω̂)
    xp = [R*x[1:3]; R*(x[4:6] - ω̂*x[1:3])]
end

function inertial_to_planet_fixed(x::AbstractVector{T},p::Planet{T}) where {T}
    #Assume axes are aligned
    inertial_to_planet_fixed(0.0,x,p)
end

function cartesian_to_vinh(x::AbstractVector{T}) where {T}
    r = norm(x[1:3])
    θ = atan(x[2], x[1]) #longitude (defined starting from x-axis in equatorial plane)
    ϕ = asin(x[3]/r) #latitude (defined up from equatoria plane)

    v = norm(x[4:6]) #velociy magnitude

    #Basis for funky LVLH frame (roughly Up, East, North)
    #See pages 2-5 and 2-6 of Hypersonic Flight Mechanics by Busemann, Vinh, and Culp
    #https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19760024112.pdf
    e1 = x[1:3]./r
    e2 = cross([0,0,1.0],e1)
    e2 = e2./norm(e2)
    e3 = cross(e1,e2)

    #Rotate velcity into LVLH frame
    v_lvlh = [e1 e2 e3]'*x[4:6]

    #flight path angle is angle between velocity and local horizontal plane
    #(plane perpendicular to r) measured positive when v rotates toward r (out)
    γ = asin(v_lvlh[1]/v)

    #heading is angle between v projected into horizontal plane and y axis of LVLH frame
    ψ = atan(v_lvlh[3],v_lvlh[2])

    x_vinh = [r, θ, ϕ, v, γ, ψ]
end

function vinh_to_cartesian(x::AbstractVector{T}) where{T}
    #unpack state vector
    r = x[1]
    θ = x[2] #longitude
    ϕ = x[3] #latitude
    v = x[4] #velocity magnitude
    γ = x[5] #flight path angle
    ψ = x[6] #heading

    #Compute position in planet-fixed frame
    sinϕ = sin(ϕ)
    cosϕ = cos(ϕ)
    sinθ = sin(θ)
    cosθ = cos(θ)
    r_cart = [r*cosϕ*cosθ, r*cosϕ*sinθ, r*sinϕ]

    #Compute velocity in a funky LVLH frame (roughly Up, East, North)
    #see pages 2-5 and 2-6 of Hypersonic Flight Mechanics by Busemann, Vinh, and Culp
    #https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19760024112.pdf
    sinγ = sin(γ)
    cosγ = cos(γ)
    sinψ = sin(ψ)
    cosψ = cos(ψ)
    v_lvlh = [v*sinγ, v*cosγ*cosψ, v*cosγ*sinψ]

    #Rotate velocity into planet-fixed frame
    e1 = r_cart./r
    e2 = cross([0,0,1.0],e1)
    e2 = e2./norm(e2)
    e3 = cross(e1,e2)
    v_cart = [e1 e2 e3]*v_lvlh

    x_cart = [r_cart; v_cart]
end
