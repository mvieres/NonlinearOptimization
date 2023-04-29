% Rosenbrock: objective and its gradient
% fv = @(x1,x2) 100 * ( x2  -  x1.^2).^2 + (1 -  x1 ).^2;                               % vectorized for fast visualization
% f  = @(x)     100 * (x(2) - x(1)^2)^2  + (1 - x(1))^2;                                % non-vectorized with vector input
% gradf = @(x)  [ -400*x(1)*x(2) + 400*x(1)^3 + 2*x(1) - 2   ;   200*(x(2)-x(1)^2) ];   % non-vectorized with vector input
% gradftol = 1d-3;

% Test function - use tol = 1.0d-6
fv    = @(x1,x2) (0.5*x1-0.2).^4 + (x2-0.6).^2        ;
f     = @(x) (0.5*x(1)-0.2)^4 + (x(2)-0.6)^2          ;  % minimum at x = [0.4 ; 0.6]
gradf = @(x) [4*0.5*(0.5*x(1)-0.2)^3  ; 2*(x(2)-0.6)] ;   
gradftol = 1d-6;

% initial point
x0{1} = [ 1.0 ; -0.5];
x0{2} = [-1.2 ;  1.0];
%x0{3} = [-10 ;  12];

% optimization settings
tol   = gradftol;  % stopping criterion || gradf || <= tol
maxit = 10000;     % maximum number of iterations

% augment function with an evaluation counter (slow, but just for the statistics)
f     = augmentWithEvalCounter(f);
gradf = augmentWithEvalCounter(gradf);

% setup parameter structures for step size strategies
armijoparams.beta    = 0.5;       % Armijo parameter: stepsize shrinking
armijoparams.gamma   = 1.0d-4;    % Armijo parameter: sufficient descent
armijoparams.minstep = 1.0d-100;  % Armijo parameter: minimum step size
wolfeparams.gamma = 1.0d-4;  % Wolfe parameter: sufficient descent
wolfeparams.eta   = 0.9;     % Wolfe parameter: curvature

% do all simulations
methods = {@armijo     , @powellwolfe};
params  = {armijoparams, wolfeparams};
k = 0;
for i = 1:length(x0)  % cycle through initial guesses x0
   for j = 1:length(methods) % cycle through step size strategies
      k = k + 1;
      f('reset');               % reset evaluation counters
      gradf('reset');           % reset evaluation counters
      X{k} = steepestdesc(f, gradf, x0{i}, tol, maxit, methods{j}, params{j});
      fn = f('count');          % retrieve evaluation counters
      gradfn = gradf('count');  % retrieve evaluation counters
      itern  = size(X{k},2)-1;  % minus 1 for x0
      desc{k} = sprintf('#%2d --> %12s:  #it = %5d', k, func2str(methods{j}), itern);
      fprintf('#%3d:  x0 = [%6.2f ; %6.2f]', k, x0{i}(1), x0{i}(2));
      fprintf('   x = [%7.4f ; %7.4f]', X{k}(1,end), X{k}(2,end));
      fprintf('   f(x) = %6.3f', f(X{k}(:,end)));
      fprintf('   %12s', func2str(methods{j}));
      fprintf('   #it = %5d   #f = %5d   #gradf = %5d', itern, fn, gradfn);
      fprintf(' \n');
   end
end


% display solution in individual figures
figbase = 333;
for k = 1:length(X)
   displaySol(figbase+k, fv, X{k}, desc{k})
end
%fprintf('Summary output of iterates:  ')
%fprintf('%g  ', cellfun(@(x) size(x,2), X))
%fprintf('\n');




