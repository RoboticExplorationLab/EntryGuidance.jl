clear

n = 3;

x4max = 103.4
rho = 1.223
c = randn(n,1)
cvx_begin 
    variable x(n)
    variable s(n)
    
    minimize(dot(c,x))
    subject to 
    norm(x) <= s
    rho*pow_pos(s,4.24) <= x4max
cvx_end