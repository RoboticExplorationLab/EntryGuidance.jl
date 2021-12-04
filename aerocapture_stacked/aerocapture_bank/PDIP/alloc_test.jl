using LinearAlgebra
using SparseArrays
using BenchmarkTools
using SuiteSparse
include(joinpath(@__DIR__,"create_MPC.jl"))
include(joinpath(@__DIR__,"lu_taylor.jl"))
struct x_cache
    c1::Vector{Float64}
    c2::Vector{Float64}
    c3::Vector{Float64}
    function x_cache(nx)
        new(zeros(nx),zeros(nx),zeros(nx))
    end
end
struct z_cache
    c1::Vector{Float64}
    c2::Vector{Float64}
    c3::Vector{Float64}
    function z_cache(nz)
        new(zeros(nz),zeros(nz),zeros(nz))
    end
end
struct y_cache
    c1::Vector{Float64}
    c2::Vector{Float64}
    c3::Vector{Float64}
    function y_cache(ny)
        new(zeros(ny),zeros(ny),zeros(ny))
    end
end
struct CACHE
    x::x_cache
    z::z_cache
    y::y_cache
    function CACHE(nx,nz,ny)
        new(x_cache(nx),z_cache(nz),y_cache(ny))
    end
end
struct IDX
    # contains the variable indexing stuff for [x;s;z;y]
    nx::Int64
    ns::Int64
    nz::Int64
    ny::Int64
    N::Int64

    x::UnitRange{Int64}
    s::UnitRange{Int64}
    z::UnitRange{Int64}
    y::UnitRange{Int64}
end
struct INIT
    LS::SparseMatrixCSC{Float64, Int64}
    # LS_F::SuiteSparse.UMFPACK.UmfpackLU{Float64, Int64}
    LS_solver::LUSparseSolver{Float64}
    rhs::Vector{Float64}
    sol::Vector{Float64}
    idx_y::UnitRange{Int64}
    function INIT(nx,nz,ny,idx,Q,G,A,idx_y)
        Ni = nx + nz + ny
        LS = spzeros(Ni, Ni)
        LS[idx.x,idx.x] = Q
        LS[idx.x,idx.s] = G'
        LS[idx.x,idx_y] = A'
        LS[idx.s,idx.x] = G
        LS[idx.s,idx.s] = -I(idx.ns)
        LS[idx_y,idx.x] = A
        new(LS,lu_sparse_solver(LS),zeros(Ni),zeros(Ni),idx_y)
    end
end
struct DELTA
    # stores all the deltas (search directions)

    # affine deltas
    x_a::Array{Float64,1}
    s_a::Array{Float64,1}
    z_a::Array{Float64,1}
    y_a::Array{Float64,1}

    # centering + correcting deltas
    x_c::Array{Float64,1}
    s_c::Array{Float64,1}
    z_c::Array{Float64,1}
    y_c::Array{Float64,1}

    # total deltas
    x::Array{Float64,1}
    s::Array{Float64,1}
    z::Array{Float64,1}
    y::Array{Float64,1}

    # constructor
    function DELTA(nx,ns,nz,ny)
        new(zeros(nx),zeros(ns),zeros(nz),zeros(ny),
            zeros(nx),zeros(ns),zeros(nz),zeros(ny),
            zeros(nx),zeros(ns),zeros(nz),zeros(ny))
    end
end

struct QP

    # problem data
    Q::SparseMatrixCSC{Float64,Int64}
    q::Array{Float64,1}
    A::SparseMatrixCSC{Float64,Int64}
    b::Array{Float64,1}
    G::SparseMatrixCSC{Float64,Int64}
    h::Array{Float64,1}

    # variables
    x::Array{Float64,1} # primal
    s::Array{Float64,1} # primal
    z::Array{Float64,1} # dual
    y::Array{Float64,1} # dual

    # KKT stuff
    KKT::SparseMatrixCSC{Float64,Int64}
    rhs_a::Array{Float64,1}
    rhs_c::Array{Float64,1}
    p_a::Array{Float64,1}
    p_c::Array{Float64,1}

    # indexing
    idx::IDX

    # deltas
    Δ::DELTA

    # cache
    cache::CACHE

    # initialization stuff
    init::INIT

    function QP(Q,q,A,b,G,h)

        # length of variables
        nx = length(q)
        ns = length(h)
        nz = length(h)
        ny = length(b)
        N = nx + ns + nz + ny

        # indexing when stacked [x;s;z;y]
        idx_x = 1:nx
        idx_s = (nx + 1) : (nx + ns)
        idx_z = (nx + ns + 1) : (nx + ns + nz)
        idx_y = (nx + ns + nz + 1) : (nx + ns + nz + ny)
        idx = IDX(nx,ns,nz,ny,N, idx_x, idx_s, idx_z, idx_y)

        # kkt stuff
        KKT = spzeros(N,N)
        rhs_a = zeros(N)
        rhs_c = zeros(N)
        p_a = zeros(N)
        p_c = zeros(N)

        # initialize variables to zero
        x = zeros(nx)
        s = zeros(ns)
        z = zeros(nz)
        y = zeros(ny)

        # deltas
        Δ = DELTA(nx,ns,nz,ny)

        #initialization
        init = INIT(nx,nz,ny,idx,Q,G,A,idx_y .- ns)
        # nx,nz,ny,idx,Q,G,A,idx_y

        new(Q,q,A,b,G,h,x,s,z,y, KKT, rhs_a, rhs_c, p_a, p_c, idx, Δ,CACHE(nx,nz,ny),init)
    end


end
function rhs_kkt_a!(qp::QP)
    idx = qp.idx
    c = qp.cache

    # qp.rhs_a[idx.x] = -(qp.A'*qp.y + qp.G'*qp.z + qp.Q*qp.x + qp.q)
    mul!(c.x.c1,qp.A',qp.y)
    mul!(c.x.c2,qp.G',qp.z)
    mul!(c.x.c3,qp.Q,qp.x)
    @. qp.rhs_a[idx.x] = -c.x.c1 - c.x.c2 - c.x.c3 - qp.q

    # qp.rhs_a[idx.s] = -(qp.z)
    @. qp.rhs_a[idx.s] = -qp.z

    # qp.rhs_a[idx.z] = -(qp.G*qp.x + qp.s - qp.h)
    mul!(c.z.c1,qp.G,qp.x)
    @. qp.rhs_a[idx.z] = -c.z.c1 - qp.s - qp.h

    # qp.rhs_a[idx.y] = -(qp.A*qp.x - qp.b)
    mul!(c.y.c1,qp.A,qp.x)
    @. qp.rhs_a[idx.y] = -c.y.c1 + qp.b
    return nothing
end
function index_sol_a!(qp::QP)
    # qp.Δ.x_a .= qp.p_a[qp.idx.x]
    qp.Δ.x_a .= view(qp.p_a,qp.idx.x)

    # qp.Δ.s_a .= qp.p_a[qp.idx.s]
    qp.Δ.s_a .= view(qp.p_a,qp.idx.s)

    # qp.Δ.z_a .= qp.p_a[qp.idx.z]
    qp.Δ.z_a .= view(qp.p_a,qp.idx.z)

    # qp.Δ.y_a .= qp.p_a[qp.idx.y]
    qp.Δ.y_a .= view(qp.p_a,qp.idx.y)
    return nothing
end
function index_sol_c!(qp::QP)
    qp.Δ.x_c .= view(qp.p_c,qp.idx.x)
    qp.Δ.s_c .= view(qp.p_c,qp.idx.s)
    qp.Δ.z_c .= view(qp.p_c,qp.idx.z)
    qp.Δ.y_c .= view(qp.p_c,qp.idx.y)
    return nothing
end
function linesearch(x,dx)
    # α = min(1.0, minimum([dx[i]<0 ? -x[i]/dx[i] : Inf for i = 1:length(x)]))
    α = 1.0
    for i = 1:length(x)
        if dx[i]<0
            α = min(α,-x[i]/dx[i])
        end
    end
    return α
end

function centering_params(qp::QP)

    μ = dot(qp.s,qp.z)/qp.idx.ns

    α = min(linesearch(qp.s,qp.Δ.s_a), linesearch(qp.z,qp.Δ.z_a))

    # σ = (dot(qp.s + α*qp.Δ.s_a, qp.z + α*qp.Δ.z_a)/dot(qp.s,qp.z))^3
    @. qp.cache.z.c1 = qp.s + α*qp.Δ.s_a
    @. qp.cache.z.c2 = qp.z + α*qp.Δ.z_a
    σ = (dot(qp.cache.z.c1,qp.cache.z.c2)/dot(qp.s,qp.z))^3

    return σ, μ
end

function rhs_kkt_c!(qp::QP, σ, μ)
    idx = qp.idx
    qp.rhs_c .= 0
    # qp.rhs_c[idx.s] = (σ*μ .- (qp.Δ.s_a .* qp.Δ.z_a)) ./ qp.s
    @. qp.rhs_c[idx.s] = (σ*μ - (qp.Δ.s_a * qp.Δ.z_a)) / qp.s
    return nothing
end
function combine_deltas!(qp::QP)
    @. qp.Δ.x = qp.Δ.x_a + qp.Δ.x_c
    @. qp.Δ.s = qp.Δ.s_a + qp.Δ.s_c
    @. qp.Δ.z = qp.Δ.z_a + qp.Δ.z_c
    @. qp.Δ.y = qp.Δ.y_a + qp.Δ.y_c
    return nothing
end
function update_vars!(qp::QP,α)
    @. qp.x += α*qp.Δ.x
    @. qp.s += α*qp.Δ.s
    @. qp.z += α*qp.Δ.z
    @. qp.y += α*qp.Δ.y
    return nothing
end
function initialize!(qp::QP)
    idx = qp.idx
    init = qp.init
    LS = qp.init.LS
    solver = qp.init.LS_solver
    # LS_F = qp.init.LS_F
    sol = qp.init.sol
    rhs = qp.init.rhs
    idx_y = qp.init.idx_y
    # @. sol = [-qp.q;qp.h;qp.b]
    @. rhs[idx.x] = -qp.q
    @. rhs[idx.s] = qp.h
    @. rhs[idx_y] = qp.b

    # factorize KKT
    factorize!(solver,LS)
    linear_solve!(solver,sol,LS,rhs,fact = true)

    @show norm(LS*sol - rhs)

    # @show norm(sol - LS\rhs)
    # ldiv!(LS,sol)
    # ldiv!(qp.init.LS_F,sol)
    # view(sol,idx.x) = -qp.q
    # Ni = idx.nx + idx.nz + idx.ny
    # @show Ni
    # idx_y = idx.y .- idx.ns
    # @show idx_y
    # @show qp.init.idx_y
    #
    # A = spzeros(Ni, Ni)
    # A[idx.x,idx.x] = qp.Q
    # LS[idx.x,idx.x] = qp.Q
    # A[idx.x,idx.s] = qp.G'
    # A[idx.x,idx_y] = qp.A'
    # A[idx.s,idx.x] = qp.G
    # A[idx.s,idx.s] = -I(idx.ns)
    # A[idx_y,idx.x] = qp.A
    #
    # init = A\[-qp.q;qp.h;qp.b]
    #
    # qp.x .= init[idx.x]
    # qp.z .= init[idx.s]
    # qp.y .= init[idx_y]
    #
    #
    # α_p = -minimum(-qp.z)
    # if α_p < 0
    #     qp.s .= -qp.z
    # else
    #     qp.s .= -qp.z .+ (1 + α_p)
    # end
    #
    # α_d = -minimum(qp.z)
    # if α_d >= 0
    #     qp.z .= qp.z .+ (1 + α_d)
    # end

    return nothing
end


function tt()
        # n = 400
        # m_eq = 250
        # m_ineq = 150

        # Q = randn(n,n);Q = Q'*Q
        # Q = sparse(Q)
        # q = randn(n)
        #
        # A = sprand(m_eq,n,0.25)
        # b = randn(m_eq)
        #
        # G = sprand(m_ineq,n,0.25)
        # h = randn(m_ineq)
        Q,q,A,b,G,h = create_MPC()
        qp = QP(Q,q,A,b,G,h)

        qp.x .= randn(length(q))
        qp.s .= randn(length(h))
        qp.z .= randn(length(h))
        qp.y .= randn(length(b))

        # @btime ff($qp)
        # @btime index_sol_c!($qp)
        # x = randn(100)
        # dx = randn(100)
        # @btime α1 = linesearch1($x,$dx)
        # @show linesearch1(x,dx)
        # @show linesearch2(x,dx)

        # @btime centering_params($qp::QP)
        # σ, μ = centering_params(qp)
        # α = 0.4

        # @btime rhs_kkt_c!($qp::QP, $σ, $μ)

        # @btime combine_deltas!($qp::QP)
        # @btime update_vars!($qp::QP,$α)
        # @btime initialize!($qp)
        initialize!(qp)


        # a = randn(10)
        # b = zeros(3)
        # idxb = 1:3
        # idx = 4:6
        # # @btime $a[$idx]
        # @show view(a,idx)
        # @btime view($a,$idx)
        # @btime $b.=view($a,$idx)
        # @btime $b = $a[$idx]
        # @btime $b .= $a[$idx]
        # @btime @. $b = $a[$idx]
        # @btime @. $b[$idxb] = $a[$idx]


        # @btime @. $b = $a[$idx]
        # @btime @. $b = view($a,$idx)


    # KKT = spzeros(50,50)
    # rhs = zeros(50)
    # idxs = 11:20
    #
    # # qp.rhs_a[idx.x] = -(qp.A'*qp.y + qp.G'*qp.z + qp.Q*qp.x + qp.q)
    # nx = 10
    # ny = 4
    # A = randn(ny,nx)
    # y = randn(ny)
    # nz = 5
    # G = randn(nz,nx)
    # z = randn(nz)
    # Q = randn(nx,nx)
    # x = randn(nx)
    #
    # q = randn(nx)
    #
    # # caches
    # c1 = zeros(nx)
    # c2 = zeros(nx)
    # c3 = zeros(nx)
    # # At = A'
    # # @btime $rhs[$idxs] = $q
    # # @btime mul!(view($rhs,$idxs) , $A', $y)
    # # @btime mul!(view($rhs,$idxs) , $A', $y)
    # # @btime mul!($c1,$A',$y)
    # # @btime mul!($c2,$G',$z)
    # # @btime mul!($c3,$Q,$x)
    # mul!(c1,A',y)
    # mul!(c2,G',z)
    # mul!(c3,Q,x)
    # @btime @. $rhs[$idxs] = $c1 + $c2 + $c3 + $q

    # @btime mul!($view($rhs,$idxs), $A', $y)

    # A = randn(5,5)
    # # B = randn(5,5)
    # B = randn(5)
    # C = randn(5)
    #
    # @btime mul!($C,$A,$B)




    return nothing
end


tt()
