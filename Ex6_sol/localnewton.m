function [X, info] = localnewton(~, gradf, hessf, x0, tol, maxit)
% X = localnewton(f, gradf, hessf, x0, tol, maxit)
%
% Locates a stationary point of function f using the local/full-step Newton method.
%
% INPUT:        f --> function handle to objective [unused!]  (accepting one input variable)
%           gradf --> function handle to objective's gradient (accepting one input variable)
%           hessf --> function handle to objective's hessian  (accepting one input variable)
%              x0 --> initial guess (starting point)
%             tol --> stopping condition tolerance:   ||gradf|| <= tol
%           maxit --> maximum number of iterations
%
% OUTPUT: 
%          X --> matrix with all iterates; X(:,k) containing the k-th iterate
%       info --> -1, if maximum number of iterations is reached
%                 0, if a stationary point (up to tol) is found
%

% setup and defaults
info     = -1;        % error status value: maximum number of iterations reached
xk       = x0;        % initial guess is zero-th iterate 
gradf_xk = inf;       % to overcome the first stopping check

% ensure initial guess x0 is a single column
x0 = reshape(x0, [], 1);

% prepare output (NOTE: Matlab uses column-major order!)
historyChunk = 1000;
dimx = length(x0);
X = zeros(dimx, historyChunk);

% store initial guess
X(:,1) = x0;
historyIdx  = 1;  % fill counter
historySize = historyChunk;

% iterate
for k = 1:maxit
   
   % check the stopping criterion 
   if norm(gradf_xk,2) <= tol
      info = 0; % set output flag
      break;
   end
   
   % compute hessian and gradient at xk
   hessf_xk = hessf(xk);  % very costly
   gradf_xk = gradf(xk);  % quite costly
   
   % calculate Newton direction
   dk = - hessf_xk \ gradf_xk;
   
   % perform the step
   xk = xk + dk;
   
   % store it in the result matrix; grow-and-copy only if necessary
   historyIdx = historyIdx + 1;
   if historyIdx > historySize
      X(:,end+1:end+historyChunk) = 0;          % increase history size
      historySize = historySize + historyChunk; % new size of history
   end
   X(:,historyIdx) = xk;
   
end

% shrink result array to appropriate size
X = X(:, 1:historyIdx);

% finito of function
return

end % end of function