%Sheet 5 main
% E4

f  = @(x) 100 * (x(2) - x(1)^2)^2 + (1-x(1))^2; 
fv = @(x1,x2) 100 * (x2 - x1.^2).^2 + (1 - x1).^2;
gradf = @(x)  [ -400*x(1)*x(2) + 400*x(1)^3 + 2*x(1) - 2   ;   200*(x(2)-x(1)^2) ];
hessf = @(x) [-400*x(2)+(1200*x(1)^2)+2, -400*x(1); -400*x(1), 200];
tol = 1e-9;
maxit = 10000;

%Calculation
x0 = [1; -0.5];
xk1 = localnewton(f,gradf, hessf, x0, tol, maxit);

x0 = [-1.2; 1];
xk2 = localnewton(f,gradf, hessf, x0, tol, maxit);

%That are suspectfully few iterates, however I did not find a bug in my
%code

x0 = [1; -0.5];
xk1_sd = steepestdesc(f,gradf,x0,tol,maxit);

x0 = [-1.2; 1];
xk2_sd = steepestdesc(f,gradf,x0,tol,maxit);


%Plots
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
title('Local Newton Method')
plot3(xk1(1,:), xk1(2,:), 0*xk1(1,:), 'k.-', 'LineWidth', 1.5, 'MarkerSize', 8.0);  
plot3(xk2(1,:), xk2(2,:), 0*xk2(1,:), 'r.-', 'LineWidth', 1.5, 'MarkerSize', 8.0);

plot3(xk1_sd(1,:), xk1_sd(2,:), 0*xk1_sd(1,:), 'b.-', 'LineWidth', 1.5, 'MarkerSize', 8.0);  
plot3(xk2_sd(1,:), xk2_sd(2,:), 0*xk2_sd(1,:), 'y.-', 'LineWidth', 1.5, 'MarkerSize', 8.0);
legend('', '','Iterate1 of Newton', 'Iterate2 of Newton', 'Iterate1 of SteepestDesc', 'Iterate2 of SteepestDesc')

