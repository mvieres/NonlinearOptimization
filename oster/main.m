%Easter sheet


load LVdata.mat
t = LVdata.t;
eta = LVdata.eta;

%% (a)
% Visulization
plot(t,eta(1,:),'blue')
hold on
plot(t,eta(2,:),'red')
hold off
legend('Prey','Predator')
title('Predator Prey')
%% (b)
% Paramteres
p_1 = [3, 2, 1, 1];
p_2 = [0.1, 0.1, 0.1, 0.1];
x0 = [6; 5];
%ode_sol = LVrhs(t,x0,p_1);

tspan = [t(1) t(end)];
[t_new, ode_sol] = ode45(@(y) odefcn(y,p_1), tspan, x0);
