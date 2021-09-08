
function rollout(model,x0,U_in,dt)
    N = length(U_in)+1
    X = [zeros(7) for i = 1:N]
    X[1] = copy(x0)
    for i = 1:(length(X)-1)
        X[i+1] = rk4(model,X[i],U_in[i],dt)

        # here we check for on the up and up
    end
    # if altitude(model,X[N])<125
    #     error("didn't make it back up")
    # end
    return X
end
