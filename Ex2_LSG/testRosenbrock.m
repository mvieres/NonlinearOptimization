% objective and its gradient
f  = @(x) 100 * (x(2) - x(1)^2)^2 + (1-x(1))^2;                                       % non-vectorized with vector input
fv = @(x1,x2) 100 * (x2 - x1.^2).^2 + (1 - x1).^2;                                    % vectorized for visualization
gradf = @(x)  [ -400*x(1)*x(2) + 400*x(1)^3 + 2*x(1) - 2   ;   200*(x(2)-x(1)^2) ];   % non-vectorized with vector input

% initial point
x0_1 = [ 1.0 ; -0.5];
x0_2 = [-1.2 ;  1.0];

% optimization settings
tol   = 1.0d-3;    % stopping criterion || gradf || <= tol
maxit = 10000;     % maximum number of iterations
gamma = 1.0d-3;    % Armijo parameter: sufficient descent
beta  = 0.5;       % Armijo parameter: stepsize shrinking


% start the optimizer
[sol1, info1, history1] = steepestdesc(f, gradf, x0_1, tol, maxit, gamma, beta);
[sol2, info2, history2] = steepestdesc(f, gradf, x0_2, tol, maxit, gamma, beta);

% visualization
figure(333); clf; hold on                    % prepare figure window
[x1, x2]  = meshgrid(-3:0.01:3, -2:0.01:4);  % get the evaluation grid
fvals     = fv(x1, x2);                      % evaluate the vectorized function
hsurf     = surf(x1, x2, fvals);             % surface plot
contours  = [0.1 1 5 10  100   1000];        % contour lines at specified f values
[~,hcont] = contour(x1, x2, fvals, contours);% contour plot
hsurf.FaceAlpha = 0.5;                       % transparency of surface plot faces
hsurf.EdgeAlpha = 0.5;                       % transparency of surface plot edges 
hcont.LineWidth = 1.5;                       % line width of contour plot
xlabel('x1');                                % set label for x-axis
ylabel('x2');                                % set label for y-axis


% add the optimization paths
xk1 = [history1.xk];  fk1 = [history1.fk];  sigmak1 = [history1.sigmak];  dk1 = [history1.dk]; % extract history
xk2 = [history2.xk];  fk2 = [history2.fk];  sigmak2 = [history2.sigmak];  dk2 = [history2.dk]; % extract history
plot3(xk1(1,:), xk1(2,:), 0*fk1, 'k.-', 'LineWidth', 1.5, 'MarkerSize', 8.0);  % plot on the z=0-plane
plot3(xk2(1,:), xk2(2,:), 0*fk2, 'r.-', 'LineWidth', 1.5, 'MarkerSize', 8.0);  % plot on the z=0-plane


% show the histories
fprintf('Result for initial guess x0 = (%f, %f)\n', x0_1(1), x0_1(2));
displayHistory(history1, 13, 6, 1, 10);   fprintf('\n\n\n');
fprintf('Result for initial guess x0 = (%f, %f)\n', x0_2(1), x0_2(2));
displayHistory(history2, 13, 6, 1, 10);   fprintf('\n\n\n');


