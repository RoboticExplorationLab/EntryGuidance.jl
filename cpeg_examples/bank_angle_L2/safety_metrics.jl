using LinearAlgebra
using EntryGuidance
function plot(x,y)
    mat"
    figure
    hold on
    plot($x,$y)
    hold off
    "
end
function plot(y)
    mat"
    figure
    hold on
    plot($y)
    hold off
    "
end
function myρ(x)
    r = x[1:3]*1000       # m
    v = x[4:6]*1000/3600  # m/s

    ρ0 = 3.032e-2*2 # kg/m^3
    H = 8.757     # km

    Rm = 3389.5 # km

    z = norm(r)/1000 - Rm # altitude km

    ρ = ρ0*exp(-z/H)
    return ρ
end

function get_pressure(model::EntryVehicle,x)
    # input is the km, km/s units
    # output is the scaled kPa units
    ρ = atmospheric_density(x[1:3], model.evmodel)*1e-9
    PR = 0.80527*ρ^(1.0036)*norm(x[4:6]*1000/3600)^2.0251
    PR /= 1000
    return PR
end
function pressure_constraint(model::EntryVehicle,x)
    # this is a less than 0 scenario
    # PR max = 20
    return get_pressure(model::EntryVehicle,x) - 18 #≦ 0
end
function pressure_gradient(model::EntryVehicle,x)
    return ForwardDiff.gradient(_x -> pressure_constraint(model,_x),x)
end
# A = ForwardDiff.gradient(_x -> pressure_constraint(M,_x),Xsim2[30])
function get_heating(model::EntryVehicle, x)

    r = x[1:3]

    # # convert v to m/s
    # v = x[4:6]*1000/3600
    # # v = x[4:6]*1/3600
    #
    # #atmospheric density # units of kg / m^3
    # ρ = atmospheric_density(r, model.evmodel)*1e-9
    #
    # # heatrate W/m^2
    # HR = 8.43e-13*ρ^(0.82958)*norm(v)^4.512

    # convert v to m/s
    v = x[4:6]*1000/3600

    #atmospheric density # units of kg / m^3
    # ρ = atmospheric_density(r, model.evmodel)*1e-9
    ρ = myρ(x)
    # ρ*= 1e6

    # heatrate W/m^2
    HR = 8.43e-13*ρ^(0.82958)*norm(v)^4.512
    # rn = 0.6 # meters
    # rn = 1
    # κ = 1.9027e-4
    #
    # HR = κ*sqrt(ρ/rn)*norm(v)^3 # per m^3
    # HR *= 1e-6

    # pressure
    # PR = 0.80527*ρ^(1.0036)*norm(v)^2.0251
    PR = get_pressure(model,x)

    # Shear stress
    # SS = 3.3e-6*ρ^(0.75356)*norm(v)^2.7409

    return HR, PR, ρ
end

M = EntryVehicle(CartesianMSLModel(),1e4)
alt_hist = zeros(length(Xsim2))
ρ_hist = zeros(length(Xsim2))
ρ_hist2 = zeros(length(Xsim2))
HR_hist = zeros(length(Xsim2))
PA_hist = zeros(length(Xsim2))
for i = 1:length(Xsim2)
    alt_hist[i] = altitude(M,Xsim2[i])
    # ρ_hist[i] = atmospheric_density(Xsim2[i][1:3],M.evmodel)*1e-9
    HR_hist[i], PA_hist[i], ρ_hist[i] = get_heating(M,Xsim2[i])
    ρ_hist2[i] = myρ(Xsim2[i])
end

accel = zeros(length(Xsim2)-1)
for i = 1:(length(Xsim2) - 1)
    v0 = Xsim2[i][4:6]*1000/3600
    v1 = Xsim2[i+1][4:6]*1000/3600
    dt = 2.0
    a = (v1-v0)/2
    # accel[i] = norm(a)/9.8
    accel[i] = norm(v1 - v0)/(2*9.8)
end
mat"
figure
hold on
plot($ρ_hist,$alt_hist)
plot($ρ_hist2,$alt_hist,'r')
hold off
set(gca, 'XScale', 'log')
"

mat"
figure
hold on
plot($HR_hist)
hold off
"

# mat"
# figure
# hold on
# plot($PA_hist)
# plot([1;length($PA_hist)],[15.4e3;15.4e3],'r')
# hold off
# "

plot(accel)

plot(PA_hist)
# plot(ρ_hist,alt_hist)
