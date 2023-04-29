function [X, info] = bfgsinverse(f, gradf, B0, x0, tol, maxit, gamma, eta, varargin)
% [X, info] = bfgsinverse(f, gradf, B0, x0, tol, maxit, gamma, eta)
% [X, info] = bfgsinverse(f, gradf, B0, x0, tol, maxit, gamma, eta, output)
%
% Inverse BFGS algorithm with Powell-Wolfe step size strategy for globalization.
%
% INPUT:        f --> function handle to objective             (accepting one input variable)
%           gradf --> function handle to objective's gradient  (accepting one input variable)
%              B0 --> initial inverse Hessian approximation
%              x0 --> initial guess (starting point)
%             tol --> stopping condition tolerance:   ||gradf|| <= tol
%           maxit --> maximum number of iterations
%           gamma --> Powell-Wolfe parameter (sufficient descent)
%             eta --> Powell-Wolfe parameter (curvature condition)
%          output --> flag: display informative output in every iteration [default: true]
%
% OUTPUT: 
%          X --> matrix with all iterates; X(:,k) containing the k-th iterate
%       info --> -1, if maximum number of iterations is reached
%                 0, if a stationary point (up to tol) is found
%

% if initial Hessian approximation is scalar, lift it to diagonal matrix
if isscalar(B0)
   B0 = B0 * eye(length(x0));
end

% setup and defaults
info     = -1;        % error status value: maximum number of iterations reached
xk       = x0;        % initial guess is zero-th iterate 
Bk       = B0;        % initial inverse Hessian approximation
output   = true;      % output status information on each iteration

% process optional input
if (nargin >= 8); output = varargin{1}; end

% ensure initial guess x0 is a single column
x0 = reshape(x0, [], 1);

% compute first gradient at xk
gradfxk = gradf(xk);

% check the stopping criterion (lucky guess?)
if norm(gradfxk,2) <= tol
   info = 0; % set output flag
   X = x0;   % "history"
   return
end

% prepare output and history arrays
historyChunk = 1000;
historySize = historyChunk;
X = zeros(length(x0), historyChunk);  % output matrix

% store initial guess
X(:,1) = x0;

% initial output
if output
   outputIteration(1, X, NaN, gradfxk, f); 
end

% iterate
for k = 1:maxit

   % search direction:
   dk = - Bk*gradfxk;
   
   % globalization by Powell Wolfe
   sigmak = powellwolfe(f, gradf, xk, dk, gamma, eta);
   
   % some helpers
   step = sigmak * dk;  % for reuse lateron
   gradfprev = gradfxk;
   
   % compute next iterate and new gradient
   xk = xk + step;
   gradfxk = gradf(xk);
   
   % store iterare in the result matrix; grow-and-copy only if necessary
   if k >= historySize
      historySize = historySize + historyChunk; % new size of history
      X(:,end+1:end+historyChunk) = 0;          % increase history size
   end
   X(:,k+1) = xk;

   % check termination criterion
   if (norm(gradfxk,2) <= tol)
      info = 0; % convergence marker
      break
   end
   
   % update inverse Hessian approximation
   sk = step;
   yk = gradfxk - gradfprev;
   Bk = updateDFP(Bk, yk, sk);

   % informative output
   if output
      outputIteration(k+1, X, sigmak, gradfxk, f); 
   end
   
   % save(sprintf('BFGS-Bk-%03d.txt',k), 'Bk', '-ascii','-tabs'); % TESTING
   
end

% shrink result array to appropriate size
X = X(:, 1:k+1);

end % end of function



