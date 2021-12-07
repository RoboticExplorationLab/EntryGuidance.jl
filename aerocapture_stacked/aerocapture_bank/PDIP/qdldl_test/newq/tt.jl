using QDLDL, LinearAlgebra, SparseArrays, SuiteSparse
using Test
function random_psd(n)

    A = sprandn(n,n,0.2)
    A = A+A';
    A = A + Diagonal(((sum(abs.(A),dims=1)[:]))) #make diagonally dominant

end

function build_kkt(Q,G)
  n = size(Q,1)
  m = size(G,1)

  A = [Q G';G zeros(m,m)] + 1e-2*Diagonal([ones(n);-ones(m)])
  A = sparse(A)

  b = randn(n + m)

  return A,b
end

function tt()
  # random KKT system
  m = 20
  n = 30
  Q = random_psd(n)
  G = sprandn(m, n, 0.3)
  A1,b1 = build_kkt(Q, G)

  # create factorisation
  F = qdldl(A1)

  # compare qdldl solve to \
  @test norm(F\b1 - A1\b1) < 1e-10

  # update kkt with new data in the same sparsity pattern
  Q2 = copy(Q)
  Q2.nzval .= randn(nnz(Q2))
  Q2 = Q2' + Q2
  G2 = copy(G)
  G2.nzval .= randn(nnz(G2))
  A2,b2 = build_kkt(Q2, G2)

  # update factorisation in-place and test solution
  update_A!(F,A2)
  @test norm(F\b2 - A2\b2) < 1e-10

  return nothing
end

tt()






# end
