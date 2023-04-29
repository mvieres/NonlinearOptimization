function [x, info, history] = steepestdesc(f, gradf, x0, tol, maxit, gamma, beta)
% x = steepestdesc(f, gradf, x0, tol, maxit)
% x = steepestdesc(f, gradf, x0, tol, maxit, gamma)
% x = steepestdesc(f, gradf, x0, tol, maxit, gamma, beta)
% [x, info] = steepestdesc(...)
% [x, info, history] = steepestdesc(...)
%
% Locates a stationary point of function f using the method of steepest descent
% with an Armijo step size strategy
%
% INPUT:   f --> function handle to objective (accepting one input variable)
%      gradf --> function handle to objective's gradient (accepting one input variable)
%         x0 --> initial guess (starting point)
%        tol --> stopping condition tolerance:   ||gradf|| <= tol
%      maxit --> maximum number of iterations
%      gamma --> Armijo-parameter for descent sufficiency  [default: 1.0d-3]
%       beta --> Armijo-parameter for step size shrinking  [default: 0.5   ]
%
% OUTPUT: 
%          x --> stationary point of f or 
%                or defined minimum step size if no sufficient descent could be found.
%       info --> -1, if maximum number of iterations is reached
%                 0, if a stationary point (up to tol) is found
%    history --> structure with following fields:
%                  .xk --> contains the xk
%                  .dk --> contains the search directions dk
%                  .sigmak --> contains the step sizes sigmak
%                  .fk --> contains the function values f(xk)
%                  .gradfk --> contains the gradients gradf(xk)

% setup and defaults
storehistory = false; % do not store the iterates
info     = -1;        % error status value: maximum number of iterations reached
xk       = x0;        % initial guess is zero-th iterate 
gradf_xk = inf;       % to overcome the first stopping check

% process optional args
if (nargin < 6), gamma = 1.0d-3; end   % Armijo parameter: sufficient descent
if (nargin < 7), beta  = 0.5;    end   % Armijo parameter: stepsize shrinking
if (nargout >= 3), storehistory = true; end   % store history only if requested

% some preparations for the first history entry
if storehistory
   fk = f(x0);
end

% iterate
for k = 1:maxit
   
   % check the stopping criterion 
   if norm(gradf_xk,2) <= tol
      info = 0; % set output flag
      break;
   end
   
   % compute gradient at xk
   gradf_xk = gradf(xk);
   
   % steepest descent search direction: negative gradient
   dk = - gradf_xk;
   % dk = dk / 2^k;   % TEST: vanishing directions
   
   % compute Armijo stepsize
   [sigmak, xnew, fnew] = armijo(f, gradf_xk, xk, dk, gamma, beta);
   
   % store the intermediates if requested
   if storehistory
      history(k) = historyentry(xk, dk, sigmak, fk, gradf_xk, xnew, fnew);
      fk = fnew;  % store the new function value (as it is f(xk) in next iteration)
      % note that history grows with every iteration, which is not optimal w.r.t. performance.
      % however, one would not store the history at all, if performance is the matter.
   end
   
   % perform the step
   % xk = xk + sigmak * dk;    % --> already done by armijo
   xk = xnew;
   
end

% set output variable
x = xk;

% finito of function
return

% some helpers
   function h = historyentry(xk, dk, sigmak, fk, gradfk, xnew, fnew)
      h = struct('xk'    , xk,     ...
                 'dk'    , dk,     ...
                 'sigmak', sigmak, ...
                 'fk'    , fk,     ...
                 'gradfk', gradfk, ...
                 'xnew'  , xnew,   ...
                 'fnew'  , fnew );
   end


% end of outer function
end