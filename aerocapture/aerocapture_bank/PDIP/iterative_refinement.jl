using LinearAlgebra


function iterative_ref(factor,A,b)
    r = copy(b)
    x = zeros(length(b))
    for i = 1:10
        x += factor\r
        r = b - A*x
        # @show norm(r)
        if norm(r)<1e-8
            break
        end
        if i == 10
            @infiltrate
            error()
        end
    end
    return x
end

#
# let
#
#     n = 10
#     A = randn(n,n)
#     b = randn(n)
#     factor = factorize(A + 1e-6*I)
#
#     x1 = A\b
#     x2 = iterative_ref(factor,A,b)
#     @show norm(x1-x2)
# end
