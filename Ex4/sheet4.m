% Sheet 4 Exercise 4

f  = @(x) 100 * (x(2) - x(1)^2)^2 + (1-x(1))^2; 
fv = @(x1,x2) 100 * (x2 - x1.^2).^2 + (1 - x1).^2;
gradf = @(x)  [ -400*x(1)*x(2) + 400*x(1)^3 + 2*x(1) - 2   ;   200*(x(2)-x(1)^2) ];

epsilon = 1e-3;
maxit = 10000;

x0 = [1; -0.5];
xk1 = steepestdesc(f,gradf,x0,epsilon,maxit);

x0 = [-1.2; 1];
xk2 = steepestdesc(f,gradf,x0,epsilon,maxit);

figure(333); clf; hold on 
[x1, x2]  = meshgrid(-3:0.01:3, -2:0.01:4);  % get the evaluation grid
fvals     = fv(x1, x2);                      % evaluate the vectorized function
hsurf     = surf(x1, x2, fvals);             % surface plot
contours  = [0.1 1 5 10  100   1000];
[~,hcont] = contour(x1, x2, fvals, contours);
hsurf.FaceAlpha = 0.5;                       % transparency of surface plot faces
hsurf.EdgeAlpha = 0.5;                       % transparency of surface plot edges 
hcont.LineWidth = 1.5;                       % line width of contour plot
xlabel('x1');                                % set label for x-axis
ylabel('x2');   
title('Powell Wolfe Line search')
plot3(xk1(1,:), xk1(2,:), 0*xk1(1,:), 'k.-', 'LineWidth', 1.5, 'MarkerSize', 8.0);  
plot3(xk2(1,:), xk2(2,:), 0*xk2(1,:), 'r.-', 'LineWidth', 1.5, 'MarkerSize', 8.0);  
