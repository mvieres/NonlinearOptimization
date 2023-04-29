f     = @(x) (2/3)*x.^3 + (1/2)*x.^2;
gradf = @(x) 2*x.^2 + x;

% initial points
x0_1 = 1;
x0_2 = 0.5;
x0_3 = 0.1;


% optimization settings
tol   = 1.0d-6;    % stopping criterion || gradf || <= tol
maxit = 10000;     % maximum number of iterations
gamma = 1.0d-3;    % Armijo parameter: sufficient descent
beta  = 0.5;       % Armijo parameter: stepsize shrinking


% start the optimizer
[sol1, info1, history1] = steepestdesc(f, gradf, x0_1, tol, maxit, gamma, beta);
[sol2, info2, history2] = steepestdesc(f, gradf, x0_2, tol, maxit, gamma, beta);
[sol3, info3, history3] = steepestdesc(f, gradf, x0_3, tol, maxit, gamma, beta);


% show the histories
fprintf('Result for initial guess x0 = %f:\n', x0_1);
displayHistory(history1, 24, 6, 1, 5);   fprintf('\n\n\n');
fprintf('Result for initial guess x0 = %f:\n', x0_2);
displayHistory(history2, 24, 6, 1, 5);   fprintf('\n\n\n');
fprintf('Result for initial guess x0 = %f:\n', x0_3);
displayHistory(history3, 24, 6, 1, 5);   fprintf('\n\n\n');
