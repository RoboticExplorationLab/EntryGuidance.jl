clear

theta = deg2rad(60);
Q = [cos(theta) -sin(theta) 0 ;
     sin(theta)  cos(theta) 0 ;
     0           0          1];
 
 
figure
hold on 
N = quiver3(zeros(3,1),zeros(3,1),zeros(3,1),[1 0 0]',[0 1 0]',[0 0 1]')
B = quiver3(zeros(3,1),zeros(3,1),zeros(3,1),Q*[1 0 0]',Q*[0 1 0]',[0 0 1]')
hold off 

 
 

