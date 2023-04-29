function [t_new, ode_sol] = LVrhs(t,x0,p)
% Syntax: ode_sol = LVrhs(t,x0,p)
%
% Input:
%       t --> Timepoint (vector)
%       x0 --> Starting Point
%       p --> parameter vector (4 entries)
%
% Output:
%       ode_sol --> Solution of Lotka Volterra

% Note: Ich hab den header dx = LVrhs(t,x,p) nicht verstanden. Vor allem
% ist x doch beim bestimmen der Lsg einer ODE kein Inputargument. Das will
% man ja gerade bekommen. Habe die Aufgabe aber auch 1h vor Abgabe gemacht,
% also habe ich wahrscheinlich was uebersehen.

% Solving the Lotka Volterra ODE
ode_f = @(t,x) [p(1)*x(t,1) - p(3)*x(t,1)*x(t,2); ...
    -p(2)*x(t,2) + p(4)*x(t,1)*x(t,2)];



tspan = [t(1) t(end)];
[t_new, ode_sol] = ode45(@(x) odefcn(y,p), tspan, x0);

end