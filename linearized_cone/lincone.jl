using LinearAlgebra, MATLAB




function Rz(θ)
    return [cos(θ) sin(θ); -sin(θ) cos(θ)]
end

function runit()

N = 100000
x = NaN*zeros(N)
y = NaN*zeros(N)
# θ = LinRange(0,pi/4,3
# θ = 0:((pi/2)/2):(pi/2)
# θ = deg2rad.(0:10:30)
θ = deg2rad.([0,45])
# @infiltrate
# error()
# nk = 8
# θ = [ (nk/2) ]
lo = -.3
hi = .3

A = zeros(length(θ)*2,2)
for k = 1:length(θ)
    A[ (k-1)*2 .+ (1:2),:] = Rz(θ[k])
end
b_lo = lo*ones(length(θ)*2)
b_hi = hi*ones(length(θ)*2)

for i = 1:N
    # @infiltrate
    # error()
    x[i] = randn()
    # @show length(x)
    y[i] = randn()
    z = A*[x[i];y[i]]

    if (sum( z .> b_hi) + sum( z .< b_lo) != 0 )
        x[i] = NaN
        y[i] = NaN
    end
    # @infiltrate
    # error()
    # x[i] = clamp(b[1],lo,hi)
    # y[i] = clamp(b[2],lo,hi)

    # if x[i]>hi
    #     error("1")
    # end
    # if x[i]<lo
    #     error("2")
    # end
    # if y[i]>hi
    #     error("3")
    # end
    # if y[i]<lo
    #     error("4")
    # end

end

mat"
figure
hold on
plot($x,$y,'*')
hold off
axis equal
"
end

runit()
