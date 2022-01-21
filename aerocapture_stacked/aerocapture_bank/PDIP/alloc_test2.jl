using LinearAlgebra
using BenchmarkTools



# function foo(x,y)
#     return dot(x,y)

function tt()


    x = randn(4)
    y = randn(4)

    @btime c = dot($x,$y)
    @btime d = norm($x)

    return nothing
end

tt()
