using LinearAlgebra, ForwardDiff, MATLAB, Infiltrator, StaticArrays
using EntryGuidance, Attitude
const EG = EntryGuidance
include(joinpath(dirname(@__DIR__),"dynamics.jl"))
include(joinpath(dirname(@__DIR__),"rollout_stuff.jl"))
# include(joinpath(dirname(@__DIR__),"mpc.jl"))
include(joinpath(dirname(@__DIR__),"post_process.jl"))

function footprint(θ)

    model = EntryVehicle(CartesianMSLModel(),1e4)
    #Initial conditions for MSL
    Rm = model.evmodel.planet.R
    r0 = [Rm+125.0, 0.0, 0.0] #Atmospheric interface at 125 km altitude
    V0 = 5.845*3600 #Mars-relative velocity at interface of 5.845 km/sec
    γ0 = -15.474*(pi/180.0) #Flight path angle at interface
    v0 = V0*[sin(γ0), cos(γ0), 0.0]
    x0 = [r0;v0]
    x0 = [3443.300786841311, 270.4345771068569, 0.0, -6051.64651501579, 20222.23824790719, 0.0]
    #Final conditions for MSL 631.979 km down range and 7.869 km cross-range from entry point
    Rf = Rm+10.0 #Parachute opens at 10 km altitude
    rf = Rf*[cos(7.869/Rf)*cos(631.979/Rf); cos(7.869/Rf)*sin(631.979/Rf); sin(7.869/Rf)]
    vf = zeros(3)
    xf = [rf;vf]

    # first rollout
    dt = 2/3600
    N = 1000
    X = NaN*[@SArray zeros(6) for i = 1:N]
    U = [@SArray zeros(2) for i = 1:N-1]

    X[1] = deepcopy(x0)
    end_idx = NaN
    for i = 1:(N-1)
        U[i] = getmaxL(model,X[i])*[sin(θ);cos(θ)]
        # U[i] = getmaxL(model,X[i])*[.01*randn();0.41]
        X[i+1] = rk4(model,X[i],U[i],dt)
        if altitude(model,X[i+1])<10
            @info "under altitude on first rollout"
            end_idx = i + 1
            break
        end
    end
    # error()
    X = X[1:end_idx]
    U = U[1:end_idx]
    # return nothing

    # P = backpass(model,X,U,xf,dt)
    # @infiltrate
    # error()

    alt, dr, cr = postprocess(model,X,x0)
    return alt, dr, cr
end


function main()

    θ_vec = LinRange(0,2*pi,100)
    T = length(θ_vec)
    althist = [zeros(2) for i = 1:T]
    drhist = [zeros(2) for i = 1:T]
    crhist = [zeros(2) for i = 1:T]
    # xf = zeros(6)
    for i = 1:T
        althist[i], drhist[i], crhist[i] = footprint(θ_vec[i])
    end

    mat"
    figure
    hold on
    rgb1 = [29 38 113]/255;
    rgb2 = [195 55 100]/255;
    drgb = rgb2-rgb1;
    for i = 1:length($drhist)
        px = $drhist{i};
        py = $crhist{i};
        alt = $althist{i};
        %plot(px,py,'Color',rgb1 + drgb*(i-1)/length($drhist),'linewidth',3)
        %plot3(px,py,alt,'b')%,'Color',rgb1 + drgb*(i-1)/length($drhist),'linewidth',3)
        plot(px,py,'b')
    end
    xlabel('Downrange (km)')
    ylabel('Crossrange (km)')
    legend('test')
    hold off
    addpath('/Users/kevintracy/devel/Low-Thrust-TrajOpt/matlab2tikz/src')
    matlab2tikz('test.tex')
    %saveas(gcf,'range.png')
    "
end
main()
