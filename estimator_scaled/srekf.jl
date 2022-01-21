using LinearAlgebra

# struct KF_SYS
#     ΓQ::UpperTriangular{Float64, Matrix{Float64}}
#     ΓR::UpperTriangular{Float64, Matrix{Float64}}
# end
function measurement(model,x)
    return copy(x[1:7])
end
function chol(A)
    # returns upper triangular Cholesky factorization of matrix A
    return cholesky(Symmetric(A)).U
end
function qrᵣ(A)
    # QR decomposition of A where only the upper triangular R is returned
    return qr(A).R
end

function sqrkalman_filter(model, μ,F,u,y,kf_sys)

    μ̄, F̄ = sqrkf_predict(model, μ,F,u,kf_sys)

    z, L = sqrkf_innovate(model, μ̄,F̄,y,kf_sys)

    μ₊, F₊ = sqrkf_update(model, μ̄,F̄,z,L,kf_sys)

    return μ₊, F₊
end

function sqrkf_predict(model, μ,F₋,u,kf_sys)

    # get probem data from kf_sys
    # A, B, ΓQ = kf_sys.A, kf_sys.B, kf_sys.ΓQ
    ΓQ = kf_sys.ΓQ
    dt = kf_sys.dt

    A = ForwardDiff.jacobian(_x -> rk4(model,_x,u,dt),μ)

    # predict one step
    # μ̄ = A*μ + B*u
    μ̄ = rk4(model,μ,u,dt)
    F̄ = qrᵣ([F₋*A';ΓQ])

    return μ̄, F̄
end

function sqrkf_innovate(model, μ̄,F̄,y,kf_sys)

    # get probem data from kf_sys
    # A, C, ΓR = kf_sys.A, kf_sys.C, kf_sys.ΓR
    ΓR = kf_sys.ΓR
    C = ForwardDiff.jacobian(_x -> measurement(model,_x),μ̄)


    # innovation
    # z = y - C*μ̄
    z = y - measurement(model,μ̄)
    G = qrᵣ([F̄*C';ΓR])

    # kalman gain
    L = ((F̄'*F̄*C')/G)/(G')

    return z, L
end
function sqrkf_update(model, μ̄,F̄,z,L,kf_sys)

    # problem data
    # ΓR, C = kf_sys.ΓR, kf_sys.C
    ΓR = kf_sys.ΓR
    C = ForwardDiff.jacobian(_x -> measurement(model,_x),μ̄)

    # update (Joseph form for Σ)
    μ₊= μ̄ + L*z
    F₊= qrᵣ([F̄*(I - L*C)';ΓR*L'])

    return μ₊, F₊
end
