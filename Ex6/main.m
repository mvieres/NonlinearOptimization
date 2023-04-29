% Sheet 6 main
diary('output_console.txt');
% Functions
f = @(x) log(exp(x(1)) + exp(-x(1))) + x(2)^2;
fv = @(x1,x2) log(exp(x1) + exp(-x2)) + x2.^2;
gradf = @(x)[ (exp(x(1))-exp(-x(1)))*(exp(x(1)) + exp(-x(1)))^-1; 2*x(2)];
hessf = @(x)[((exp(x(1))+exp(-x(1)))^2 - (exp(x(1))+exp(-x(1)))* ...
    (exp(x(1))-exp(-x(1))))/((exp(x(1))+exp(-x(1)))^2), 0; 0, 2];

%gradf = @(x)  [  (exp(2*x(1))-1)/(exp(2*x(1))+1) ;   2*x(2) ];   
%hessf = @(x) [[4*exp(2*x(1))/(exp(2*x(1))+1)^2, 0];[0, 2]];

% Parameters
epsilon = 1e-6;
maxit = 100;
beta = 0.5;
gamma = 1e-4;
a1 = 1e-6;
a2 = a1;
p = 0.1;


% Computations
x_a = [1 ; -0.5];
xk1 = globalnewton(f, gradf, hessf, x_a, ...
    epsilon, maxit, beta, gamma, a1, a2, p);
xk1_l = localnewton(gradf, hessf, x_a, epsilon, maxit);

x_b = [-1.2; 1];
xk2 = globalnewton(f, gradf, hessf, x_b, ...
    epsilon, maxit, beta, gamma, a1, a2, p);
xk2_l = localnewton(gradf, hessf, x_b, epsilon, maxit);

x_c = [10; 12];
xk3 = globalnewton(f, gradf, hessf, x_c, ...
    epsilon, maxit, beta, gamma, a1, a2, p);
%xk3_l = localnewton(gradf, hessf, x_c, epsilon, maxit);

x_d = [33; -62];
xk4 = globalnewton(f, gradf, hessf, x_d, ...
    epsilon, maxit, beta, gamma, a1, a2, p);
%xk4_l = localnewton(gradf, hessf, x_d, epsilon, maxit);

% Note: Local newton does not work for x_c and x_d since the Hessian is
% singular at some points, as one can observe when looking at GNrGNF etc.

% Plots
figure(1); clf; hold on
[x1, x2]  = meshgrid(-10:0.5:35, -63:0.5:15); 
fvals     = fv(x1, x2);                      
hsurf     = surf(x1, x2, fvals);            
contours  = [0.001 0.05 0.1 0.5 1 3 5 10 20 30 100 1000];
[~,hcont] = contour(x1, x2, fvals, contours);
hsurf.FaceAlpha = 0.5;                       
hsurf.EdgeAlpha = 0.5;                        
hcont.LineWidth = 1.5;                       
xlabel('x1');                                
ylabel('x2');   
title('Globalized Newton Method')
plot3(xk1(1,:), xk1(2,:), 0*xk1(1,:), ...
    'k.-', 'LineWidth', 1.5, 'MarkerSize', 10.0); 
plot3(xk1_l(1,:), xk1_l(2,:), 0*xk1_l(1,:), ...
    '--xk', 'LineWidth', 1.5, 'MarkerSize', 10.0); 
plot3(xk2(1,:), xk2(2,:), 0*xk2(1,:), ...
    'g.-', 'LineWidth', 1.5, 'MarkerSize', 8.0); 
plot3(xk2_l(1,:), xk2_l(2,:), 0*xk2_l(1,:), ...
    '--xg', 'LineWidth', 1.5, 'MarkerSize', 10.0); 
plot3(xk3(1,:), xk3(2,:), 0*xk3(1,:), ...
    'r.-', 'LineWidth', 1.5, 'MarkerSize', 8.0); 

plot3(xk4(1,:), xk4(2,:), 0*xk4(1,:), ...
    'm.-', 'LineWidth', 1.5, 'MarkerSize', 8.0); 

legend(' ',' ','x_a Global','x_a local', 'x_b Global','x_b local',...
    'x_c Global', 'x_d Global');
hold off

diary("off")