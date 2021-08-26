
function rollout(model,x0,U_in,dt)
    N = 1000
    X = [@SVector zeros(6) for i = 1:N]
    U = [@SVector zeros(2) for i = 1:N-1]
    X[1] = copy(x0)
    end_idx = NaN
    t_vec = 0:dt:((N-1)*dt)
    for i = 1:(length(X)-1)

        if i > length(U_in)
            U[i] = deepcopy(U_in[end])
            # @info "went past planned controls"
        else
            U[i] = deepcopy(U_in[i])
        end

        # step forward in sim
        X[i+1] = rk4(model,X[i],U[i],dt)

        # check for hitting the 10,000ft mark
        if altitude(model,X[i+1]) < 10
            end_idx = i+1
            break
        end
    end
    # trim relavent
    if isnan(end_idx)
        error("didn't hit the altitude during the rollout")
    end

    X = X[1:end_idx]
    U = U[1:(end_idx-1)]
    t_vec = t_vec[1:end_idx]

    t_impact = 0
    return X, U, t_vec, t_impact
end
