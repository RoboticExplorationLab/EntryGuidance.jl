using LinearAlgebra
using SuiteSparse
using SparseArrays
using BenchmarkTools

"""
    LU (sparse) solver
"""
mutable struct LUSparseSolver{T}
    A::SparseMatrixCSC{T,Int}
    F::SuiteSparse.UMFPACK.UmfpackLU{T, Int}
    umf_ctrl::Vector{T}
    umf_info::Vector{T}
    tmp::Vector{Ptr{Cvoid}}
end

function lu_sparse_solver(A::SparseMatrixCSC{T,Int}) where T
    F = lu(copy(A))
    UMFPACK_INFO = 90
    UMFPACK_CONTROL = 20
    umf_ctrl = Vector{Float64}(undef, UMFPACK_CONTROL)
    umf_info = Vector{Float64}(undef, UMFPACK_INFO)
    tmp = Vector{Ptr{Cvoid}}(undef, 1)
    LUSparseSolver(copy(A), F, umf_ctrl, umf_info, tmp)
end

function factorize!(s::LUSparseSolver{T}, A::SparseMatrixCSC{T,Int}) where T
    _lu!(s.F, A, s.tmp, s.umf_ctrl, s.umf_info)
    # s.F = lu(A)
end

function linear_solve!(s::LUSparseSolver{T}, x::Vector{T}, A::SparseMatrixCSC{T,Int},
        b::Vector{T}; reg::T = 0.0, fact::Bool = true) where T
    fact && factorize!(s, A)
    _solve!(x, s.F, b, s.umf_ctrl, s.umf_info)
    # x .= A \ b
end

function linear_solve!(s::LUSparseSolver{T}, x::Matrix{T}, A::Matrix{T},
    b::Matrix{T}; reg::T = 0.0, fact::Bool = true) where T
    fill!(x, 0.0)
    n, m = size(x)
    r_idx = 1:n
    fact && factorize!(s, A)
    for j = 1:m
        xv = @views x[r_idx, j]
        bv = @views b[r_idx, j]
        _solve!(xv, s.F, bv, s.umf_ctrl, s.umf_info)
    end
end

function _umfpack_numeric!(U::SuiteSparse.UMFPACK.UmfpackLU{Float64,Int64}, tmp, umf_ctrl, umf_info; reuse_numeric = true)
    status = ccall((:umfpack_di_numeric, :libumfpack), Int64,
                   (Ptr{Int64}, Ptr{Int64}, Ptr{Float64}, Ptr{Cvoid}, Ptr{Cvoid},
                    Ptr{Float64}, Ptr{Float64}),
                   U.colptr, U.rowval, U.nzval, U.symbolic, tmp,
                   umf_ctrl, umf_info)
    return nothing
end

function _lu!(F::SuiteSparse.UMFPACK.UmfpackLU, S::SparseMatrixCSC{<:SuiteSparse.UMFPACK.UMFVTypes,<:SuiteSparse.UMFPACK.UMFITypes}, tmp, umf_ctrl, umf_info; check::Bool=true)
    zerobased = SparseArrays.getcolptr(S)[1] == 0
    F.m = size(S, 1)
    F.n = size(S, 2)

    F.colptr .= SparseArrays.getcolptr(S)
    F.rowval .= SparseArrays.rowvals(S)
    if !zerobased
        F.colptr .-= 1
        F.rowval .-= 1
    end

    F.nzval .= nonzeros(S)

    _umfpack_numeric!(F, tmp, umf_ctrl, umf_info, reuse_numeric = false)
end

function _solve!(x, lu::SuiteSparse.UMFPACK.UmfpackLU{T, Int}, b, umf_ctrl, umf_info) where T
    ccall((:umfpack_di_solve, :libumfpack), Int64,
        (Int64, Ptr{Int64}, Ptr{Int64}, Ptr{Float64},
         Ptr{Float64}, Ptr{Float64}, Ptr{Cvoid}, Ptr{Float64},
         Ptr{Float64}),
        0, lu.colptr, lu.rowval, lu.nzval,
        x, b, lu.numeric, umf_ctrl,
        umf_info)
    return nothing
end

function ttt()

    n = 10
    S = sprand(n,n,1.0)
    b = randn(n)
    x = randn(n)

    solver= lu_sparse_solver(S)
    @show b
    linear_solve!(solver,x,S,b)

    @show b
    # n = 100
    # LS = sprand(n,n,0.3) + I
    # rhs = randn(n)
    # sol = zeros(n)
    #
    # solver = lu_sparse_solver(LS)
    #
    # factorize!(solver,LS)
    #
    # linear_solve!(solver,sol,LS,rhs)
    #
    # @show norm(LS*sol - rhs)
    #
    # @show norm(sol)
    # @show norm(rhs)
    # @show norm(LS\rhs)

    # # @btime linear_solve!($solver, $sol, $A, $b)
    # # @btime $A\$b
    # # @show norm(sol - A\b)
    #
    #
    # A2 = 2*A
    # # solver.A = copy(A2)
    # @btime factorize!($solver,$A2)
    # factorize!(solver,A2)
    # @btime linear_solve!($solver,$sol,$A2,$b,fact = false)
    # linear_solve!(solver,sol,A2,b,fact = false)
    #
    # @show norm(sol - A2\b)
    return nothing
end

ttt()
