using Attitude, LinearAlgebra, ForwardDiff, MATLAB
const FD = ForwardDiff


function spin(x)
    J = Diagonal([1;2;3])
    return J\(-cross(x,J*x))
end
function srk4(x,dt)
    k1 = dt*spin(x)
    k2 = dt*spin(x + k1/2)
    k3 = dt*spin(x + k2/2)
    k4 = dt*spin(x + k3)
    return x + (1/6)*(k1 + 2*k2 + 2*k3 + k4)
end

function midpoint(x,dt)
    # x2 = copy(x)
    x2 = x + dt*spin(x)
    # J = factorize(FD.jacobian(_x2 -> (dt*spin(0.5*(x + _x2)) - _x2), x2))

    # chist = NaN*zeros(10)
    iter = 0
    for i = 1:10
        r = x + dt*spin(0.5*(x + x2)) - x2
        # chist[i] = norm(r)
        J = FD.jacobian(_x2 -> (dt*spin(0.5*(x + _x2)) - _x2), x2)
        x2 += -(J\r)
        if norm(r)<1e-12
            iter = i
            # @info "success"
            # @show i
            break
        end
    end
    return x2 , iter
end
function getB(x,dt)
    x2 = x + dt*spin(x)
    B = (FD.jacobian(_x2 -> (dt*spin(0.5*(x + _x2)) - _x2), x2))
    return B
end
function testB(B)
    for j = 1:4
        B .+= randn(size(B))
    end
    return nothing
end
function bmidpoint(x,dt)
    x2 = x + dt*spin(x)
    # x2 = copy(x)
    B = (FD.jacobian(_x2 -> (dt*spin(0.5*(x + _x2)) - _x2), x2))

    # chist = NaN*zeros(10)
    iter = 0
    for i = 1:10
        r = x + dt*spin(0.5*(x + x2)) - x2
        # chist[i] = norm(r)
        s = -(B\r)
        x2 += s
        y = (x + dt*spin(0.5*(x + x2)) - x2) - r
        B += ((y - B*s)*s')/dot(s,s)
        if norm(r)<1e-12
            # @info "success"
            # @show i
            iter = i
            break
        end
    end
    return x2 , iter
end
function linesearch(x,x2,dt,p)
    α = 1.0
    # nr = norm( x + dt*spin(0.5*(x + x + α*p)) - (x + α*p))
    # while
    # x2 = x + α*p
    nr = norm(x + dt*spin(0.5*(x + x2)) - x2)
    for i = 1:10
        x2_t = x2 + α*p
        nr_t = norm(x + dt*spin(0.5*(x + x2_t)) - x2_t)
        if nr_t < nr
            break
        else
            α*=0.5
        end
        if i == 10
            error("linesearch failed")
        end
    end
    return α
end

    # r = x + dt*spin(0.5*(x + x2)) - x2
    # for i = 1:10
        # x2 = x + α*p
    #     r = x + dt*spin(0.5*(x + x2)) - x2
function bmidpoint2(B,x,dt)
    x2 = x + dt*spin(x)
    # x2 = copy(x)
    # B = (FD.jacobian(_x2 -> (dt*spin(0.5*(x + _x2)) - _x2), x2))

    # chist = NaN*zeros(10)
    iter = 69
    for i = 1:200
        r = x + dt*spin(0.5*(x + x2)) - x2
        # chist[i] = norm(r)
        p = -(B\r)
        s = linesearch(x,x2,dt,p)*p
        # s = -0.5*(B\r)
        x2 += s
        y = (x + dt*spin(0.5*(x + x2)) - x2) - r
        B .+= ((y - B*s)*s')/dot(s,s)
        if norm(r)<1e-12
            # @info "success"
            # @show i
            iter = i
            break
        end
    end
    return x2 , iter
end
let

    x0 = [0.2;7;0.2]

    dt = 0.1
    N = 100
    X = [zeros(3) for i = 1:N]
    ihist = NaN*zeros(N)
    X[1] = copy(x0)

    B = getB(x0,dt)
    @show B

    for i = 1:(N-1)
        # X[i+1] = srk4(X[i],dt)
        # X[i+1], ihist[i] = midpoint(X[i],dt)
        # B1 = copy(B)
        X[i+1], ihist[i] = bmidpoint2(B,X[i],dt)
        # # testB(B)
        # B2 = copy(B)
        # @show norm(B1-B2)
        # error()
        # error()
    end

    Xm = mat_from_vec(X)

    mat"
    figure
    hold on
    plot($Xm')
    hold off
    "

    mat"
    figure
    hold on
    title('iteration count for midpoint convergence')
    histogram($ihist(1:end-1))
    hold off
    "
    @show ihist
end
