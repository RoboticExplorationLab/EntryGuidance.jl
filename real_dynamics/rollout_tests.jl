
function testroll()


x0 = [-1200,-1200,1000,130,150,-20.0]
U_in = [[-0;0] for i = 1:1000]
dt = 0.5
X, U, t_vec, t_impact = rollout(x0,U_in,dt)
xm = mat_from_vec(X)
@show length(X)
mat"
figure
hold on
plot($xm(1:3,:)')
hold off
"


mat"
figure
hold on
plot3($xm(1,:),$xm(2,:),$xm(3,:))
plot3([0],[0],[0],'r.','markersize',30)
axis equal
hold off
grid on
"


end


testroll()
