using MATLAB
using ForwardDiff
using BenchmarkTools
using StaticArrays

struct DensityParameters
    a::Float64
    b::Float64
    c::Float64
    d::Float64
    e::Float64
    f::Float64
    g::Float64
    h::Float64
    i::Float64
end


# -3.994917108845906
# -0.015077714054650
# -0.125933784380101
# -0.000288818937338
#  0.011301477449896
#  0.000040799036899
# -0.000433754775214
# -0.000000176283463
#  0.000000433007741


 function density(p::DensityParameters, h)
    h = h / 1000
    if h > 125.0
        h = 125.0
    elseif h < 0.2
        h = 0.2
    end
    num = @evalpoly(h, p.a, p.c, p.e, p.g, p.i)
    den = @evalpoly(h, 1.0, p.b, p.d, p.f, p.h)
    # @show num
    # @show den
    (num/den)
end

# let
#     params = [  -4.001833776317166
#   -0.015696377977412
#   -0.129412788692219
#    0.000199820058253
#    0.010377518080309
#    0.000043652882189
#   -0.000539682362508
#   -0.000000205086106
#    0.000000874980179]
#     p = DensityParameters(params...)
#
#     # @show density(p,20)
#     hs = 1e3*Vector(3:130)
#     #
#     lrho = [density(p,hs[i]) for i =1:length(hs)]
#     #
#     mat"
#     figure
#     hold on
#     plot($hs,$lrho)
#     hold off
#     "
#
#
# end

GM_MARS = 4.2828375816e13
J2_MARS = 1960.45e-6
R_MARS = 3386.2

function j2_1(r)
        # eci position stuff
    r = norm(r)
    x,y,z = r

    # precompute repeated stuff
    Re_r_sqr = 1.5*J2_mars*(R_MARS/r)^2
    five_z_sqr = 5*z^2/r^2

return  (-μ_mars/r^3)*[x*(1 - Re_r_sqr*(five_z_sqr - 1));
                  y*(1 - Re_r_sqr*(five_z_sqr - 1));
                  z*(1 - Re_r_sqr*(five_z_sqr - 3))]
end
# function ecef_Q_ned_mat(longitude,latitude)
#     ϕ = latitude
#     λ = longitude
#
#     ecef_Q_ned = [-sin(ϕ)*cos(λ)  -sin(λ)  -cos(ϕ)*cos(λ);
#                   -sin(ϕ)*sin(λ)   cos(λ)  -cos(ϕ)*sin(λ);
#                    cos(ϕ)          0.0     -sin(ϕ)]
#
#     return ecef_Q_ned
# end
# function j2_2(r)
