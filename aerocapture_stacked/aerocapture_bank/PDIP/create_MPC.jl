function create_MPC()
    nx = 6
    nu = 3
    N= 10
    xg = [-1;-2;-3;0;0;0]
    nz=(N*nx) + (N-1)*nu
    Ac = Array([zeros(3,3) I(3);zeros(3,6)])
    Bc = Array([zeros(3,3);I(3)])

    dt = 0.5
    H = exp(dt*[Ac Bc; zeros(3,9)])
    Ad = sparse(H[1:nx,1:nx])
    Bd = sparse(H[1:nx,nx+1:end])

    xi = [(nu+nx)*(i-1) .+ (1:nx)  for i = 1:N]

    ui = [(nu+nx)*(i-1) .+ (nx .+ (1:nu))  for i = 1:N-1]
    ci = [(nx)*(i-1) .+ (1:nx)  for i = 1:N]

    # cost function
    Qc = I(nx)
    Rc = 5*I(nu)
    Q = spzeros(nz,nz)
    q = zeros(nz)
    for i = 1:N
        if i<N
            Q[xi[i],xi[i]] = Qc
            q[xi[i]] = -Qc*xg
            Q[ui[i],ui[i]] = Rc
        else
            Q[xi[i],xi[i]] = Qc
            q[xi[i]] = -Qc*xg
        end
    end


    A = spzeros(N*nx,nz)
    for i = 1:N-1
        # xkp1 - axk - buk
        A[ci[i],xi[i+1]] = I(nx)
        A[ci[i],xi[i]] = -Ad
        A[ci[i],ui[i]] = -Bd
    end
    A[ci[end],xi[1]] = I(nx)
    b = zeros((N*nx))
    b[ci[end]] = [1;4;-8;.1;.2;-.3]

    G_up = spzeros(nz,nz)
    gi = [(nu)*(i-1) .+ (1:nu)  for i = 1:N-1]
    for i = 1:N-1
        G_up[gi[i],ui[i]] = I(nu)
    end
    G = [G_up;-G_up]
    h = 2*ones(2*nz)



    return Q,q,A,b,G,h
end
