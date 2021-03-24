

function find_impact(t,z,z_itp)

    z1 = z[end-1]

    t1 = t[end-1]
    t2 = t[end]

    tg = t1
    Δt = 1e-5
    for i = 1:5
        zg = z_itp(tg)
        dz = (z_itp(tg + Δt) - z_itp(tg))/Δt
        tg -= zg/dz

        # error checking and convergence
        if tg > t2 || tg<t1
            @error "impact time is outside acceptable range"
        end
        if abs(zg)<1e-5
            break
        end
        if i == 5
            @warn "max iters reached on impact timing"
        end
    end
    return tg
end

function rollout(x0,U_in,dt)
    N = 1000
    X = [zeros(6) for i = 1:N]
    U = [zeros(2) for i = 1:N-1]
    X[1] = copy(x0)
    end_idx = NaN
    t_vec = 0:dt:((N-1)*dt)
    for i = 1:(length(X)-1)

        if i > length(U_in)
            U[i] = deepcopy(U_in[end])
            @warn "went past planned controls"
        else
            U[i] = deepcopy(U_in[i])
            # U[i] = deepcopy(U_in[end])
        end

        # step forward in sim
        X[i+1] = rk4(X[i],U[i],dt)

        # check for hitting the ground
        if X[i+1][3] < 0
            end_idx = i+1
            break
        end
    end
    # trim relavent
    if isnan(end_idx)
        error("didn't hit the ground during the rollout")
    end

    X = X[1:end_idx]
    U = U[1:(end_idx-1)]
    t_vec = t_vec[1:end_idx]

    # find impact time
    xm = mat_from_vec(X)
    z_itp = CubicSplineInterpolation(t_vec,xm[3,:])
    t_impact = find_impact(t_vec,xm[3,:],z_itp)
    return X, U, t_vec, t_impact
end
