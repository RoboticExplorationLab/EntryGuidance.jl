using MATLAB
using ForwardDiff

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

let
    params = [  -4.001833776317166
  -0.015696377977412
  -0.129412788692219
   0.000199820058253
   0.010377518080309
   0.000043652882189
  -0.000539682362508
  -0.000000205086106
   0.000000874980179]
    p = DensityParameters(params...)

    # @show density(p,20)
    hs = 1e3*Vector(3:130)
    #
    lrho = [density(p,hs[i]) for i =1:length(hs)]
    #
    mat"
    figure
    hold on
    plot($hs,$lrho)
    hold off
    "


end
