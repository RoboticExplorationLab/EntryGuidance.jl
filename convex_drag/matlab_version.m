clear
% parameters
Cd = 2.2;
Area = 4.5;
rho = 1.22e-5;
% # g = zeros(3)
g = [0;0;-3.71];
% g = zeros(3,1);


% # discrete dynamics for double integrator
dt = 0.1;
Ad = eye(6);
Ad(1:3,4:6)=dt*eye(3);
Bd = [.5*dt^2*eye(3);eye(3)];

% # initial conditions
rk = [1000;2000;3000];
vk = [500;500;10];
N = 60;

cvx_begin
    variable d(3,N-1)
    variable Gamma(N-1)
    variable state(6,N)
    
    minimize(sum_square(vec(state(4:6,:))) + 1e6*sum(Gamma))
%     p = 0;
%     for i = 1:N
%         p = p + .5*sum_square(state(4:6,i)) + 3.71*state(3,i);
%     end
%     minimize(.5*sum_square(vec(state(4:6,end))) + 3.71*state(3,end) + 1e6*sum(Gamma))
%     minimize( p + 1e3*sum(Gamma))
    
    subject to 
        
        % initial conditions 
        state(1:3,1) == rk 
        state(4:6,1) == vk 
        
        for i = 1:N-1
            state(:,i+1) == Ad*state(:,i) + Bd*(g + d(:,i))
        end
        
        for i = 1:N-1
            norm(d(:,i)) <= Gamma(i)
            .5*rho*Cd*Area*sum_square(state(4:6,i)) <= Gamma(i)
        end
        
cvx_end
        


r_cvx = state(1:3,:);
v_cvx = state(4:6,:);
d_cvx = d;
gamma_cvx = Gamma;


%%

% # now we get the true solution
r_tru = zeros(3,N);
v_tru = zeros(3,N);
d_tru = zeros(3,N-1);
r_tru(:,1) = rk;
v_tru(:,1) = vk;
Gamma_tru = zeros(N-1);
for i = 1:N-1
    vl = v_tru(:,i);
    Gamma_tru(i) = .5*rho*Cd*Area*dot(vl,vl);
    d_tru(:,i) = -.5*rho*Cd*Area*dot(vl,vl)*(vl)/norm(vl);
    next_state = Ad*[r_tru(:,i);v_tru(:,i)] + Bd*(d_tru(:,i) + g);
    r_tru(:,i+1) = next_state(1:3);
    v_tru(:,i+1) = next_state(4:6);
end


%% 
figure
hold on
quiver3(r_cvx(1,1:end-1),r_cvx(2,1:end-1),r_cvx(3,1:end-1),d_cvx(1,:),d_cvx(2,:),d_cvx(3,:))
quiver3(r_tru(1,1:end-1),r_tru(2,1:end-1),r_tru(3,1:end-1),d_tru(1,:),d_tru(2,:),d_tru(3,:))
legend('CVX','Sim')
view(-15,5)
hold off
        


%% see how well slack is working 
% dnorm_tru = [norm(d_tru[:,i]) for i = 1:size(d_tru,2)]
% dnorm_cvx = [norm(d_cvx[:,i]) for i = 1:size(d_cvx,2)]

dnorm_tru = zeros(size(d_tru,2),1);
dnorm_cvx = zeros(size(d_cvx,2),1);
real_d_cvx = zeros(size(d_cvx,2),1);
for i = 1:size(d_tru,2)
    dnorm_tru(i) = norm(d_tru(:,i));
    dnorm_cvx(i) = norm(d_cvx(:,i));
    v = v_cvx(:,i);
    real_d_cvx(i) = .5*rho*Cd*Area*dot(v,v);
end
    
figure
hold on
title('Drag Norms from CVX')
plot(real_d_cvx(1:end-1))
plot(gamma_cvx,'o')
plot(dnorm_cvx,'--')
hold off
