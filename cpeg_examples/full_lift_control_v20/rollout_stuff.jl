
function rollout(model,x0,U_in,dt)
    N = 1000
    X = [zeros(6) for i = 1:N]
    U = [zeros(2) for i = 1:N-1]
    X[1] = copy(x0)
    end_idx = NaN
    for i = 1:(length(X)-1)
        U[i] = i > length(U_in) ? deepcopy(U_in[end]) : deepcopy(U_in[i])

        # step forward in sim
        X[i+1] = rk4(model,X[i],U[i],dt)

        # check for hitting the 10,000ft mark
        if altitude(model,X[i+1][1:3]*model.dscale) < 10
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
    return X, U
end
