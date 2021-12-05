# using QDLDL, LinearAlgebra, SparseArrays, SuiteSparse
include(joinpath(@__DIR__,"qdldl_myown.jl"))
using BenchmarkTools
using LinearAlgebra

function _permute_symmetric2(A::SparseMatrixCSC{Tv, Ti}, iperm::AbstractVector{Ti}, Pr::AbstractVector{Ti}, Pc::AbstractVector{Ti}, Pv::AbstractVector{Tv},num_entries::AbstractVector{Ti}) where {Tv <: AbstractFloat, Ti <: Integer}
    # 1. count number of entries that each column of P will have
    n = size(A, 2)
    # @show num_entries
    # num_entries = zeros(Ti, n)
    num_entries .= 0
    # @show num_entries
    # error()
    Ar = A.rowval
    Ac = A.colptr
    Av = A.nzval
    # count the number of upper-triangle entries in columns of P, keeping in mind the row permutation
    for colA = 1:n
        colP = iperm[colA]
        # loop over entries of A in column A...
        for row_idx = Ac[colA]:Ac[colA+1]-1
            rowA = Ar[row_idx]
            rowP = iperm[rowA]
            # ...and check if entry is upper triangular
            if rowA <= colA
                # determine to which column the entry belongs after permutation
                col_idx = max(rowP, colP)
                num_entries[col_idx] += one(Ti)
            end
        end
    end
    # 2. calculate permuted Pc = P.colptr from number of entries
    Pc[1] = one(Ti)
    @inbounds for k = 1:n
        Pc[k + 1] = Pc[k] + num_entries[k]

        # reuse this vector memory to keep track of free entries in rowval
        num_entries[k] = Pc[k]
    end
    # use alias
    row_starts = num_entries

    # 3. permute the row entries and position of corresponding nzval
    for colA = 1:n
        colP = iperm[colA]
        # loop over rows of A and determine where each row entry of A should be stored
        for rowA_idx = Ac[colA]:Ac[colA+1]-1
            rowA = Ar[rowA_idx]
            # check if upper triangular
            if rowA <= colA
                rowP = iperm[rowA]
                # determine column to store the entry
                col_idx = max(colP, rowP)

                # find next free location in rowval (this results in unordered columns in the rowval)
                rowP_idx = row_starts[col_idx]

                # store rowval and nzval
                Pr[rowP_idx] = min(colP, rowP)
                Pv[rowP_idx] = Av[rowA_idx]

                # increment next free location
                row_starts[col_idx] += 1
            end
        end
    end
    P = SparseMatrixCSC{Tv, Ti}(n, n, Pc, Pr, Pv)
    # order row indices within P.rowcal[P.colptr[k]:P.colptr[k+1]-1]
    return (P')'
    # return nothing
end


function tt()


    n = 5
    m = 3
    Q = randn(n,n)
    Q = Q'*Q + I
    G = randn(m,n)

    A = [Q G';G zeros(m,m)]

    ρ = 1e-2
    A += Diagonal([ρ*ones(n);-ρ*ones(m)])
    A = sparse(A)
    b = randn(n + m)
    # x1 = zeros(n + m)
    #
    # @. x1 = b
    #
    # F = qdldl(A;perm = nothing)
    #
    # # @btime factor!($F.workspace,false)
    # A2 = pi*A
    # # F.workspace.triuA .= triu(A2)
    # # factor!(F.workspace,false)
    #
    # # x2 = F\b
    # x2 = zeros(length(b))
    # @. x2 = b
    # solve!(F,x2)
    #
    # @show norm(x2 - A2\b)

    p = amd(A)
    iperm = invperm(p)
    # @show p

    x1 = zeros(n + m)

    @. x1 = b

    F = qdldl(A;perm = p)


    # Pr = zeros(Int,nnz(A))
    # Pc = zeros(nnz(A))
    # Pv = zeros(nnz(A))

    Pr = zeros(Int64, nnz(A))
    Pc = zeros(Int64, size(A, 1) + 1)
    Pv = zeros(Float64, nnz(A))

    n = size(A,2)
    num_entries = zeros(Int64, n)

    # @btime P = _permute_symmetric2($A,$iperm,$Pr,$Pc,$Pv,$num_entries)
    _permute_symmetric2(A,iperm,Pr,Pc,Pv,num_entries)

    A2 = pi*A
    F.workspace.triuA .= _permute_symmetric2(A2,iperm,Pr,Pc,Pv,num_entries)
    factor!(F.workspace,false)

    # x2 = F\b
    x2 = zeros(length(b))
    @. x2 = b
    solve!(F,x2)

    @show norm(x2 - A2\b)



end

tt()
