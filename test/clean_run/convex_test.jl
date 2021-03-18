using Convex, Mosek, MosekTools




xx = Variable(3,10)

Q = randn(3,3); Q = Q'*Q;

y = normalize(randn(3))
p = norm(cholseky(Q.U)*(I - y*y')*xx[:,1],2)

p = 0
p += quadform(xx[:,2],Q)
