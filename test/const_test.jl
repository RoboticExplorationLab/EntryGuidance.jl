using LinearAlgebra, BenchmarkTools


normal_vec = randn(3000)
const const_vec = randn(3000)

function read_var(j)
    a = 0
    for i = 1:j
        a += normal_vec[i]
    end
    a
end
@btime read_var(4)
function read_const(j)
    a = 0
    for i = 1:j
        a += const_vec[i]
    end
    a
end
@btime read_const(4)
