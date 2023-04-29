function [X, info] = globalnewton(f, gradf, hessf, x0, tol, maxit, beta, gamma, a1, a2, p)
% X = globalnewton(f, gradf, hessf, x0, tol, maxit, beta, gamma, a1, a2, p)
%
% Locates a stationary point of function f using the globalized Newton method
% (Armijo step size strategy with gradient descent fallback)
%
% INPUT:        f --> function handle to objective            (accepting one input variable)
%           gradf --> function handle to objective's gradient (accepting one input variable)
%           hessf --> function handle to objective's hessian  (accepting one input variable)
%              x0 --> initial guess (starting point)
%             tol --> stopping condition tolerance:   ||gradf|| <= tol
%           maxit --> maximum number of iterations
%            beta --> Armijo parameter (step size shrinking)
%           gamma --> Armijo parameter (sufficient descent)
%       a1, a2, p --> parameters for direction selection
%
% OUTPUT: 
%          X --> matrix with all iterates; X(:,k) containing the k-th iterate
%       info --> -1, if maximum number of iterations is reached
%                 0, if a stationary point (up to tol) is found
%


% setup and defaults
info     = -1;        % error status value: maximum number of iterations reached
xk       = x0;        % initial guess is zero-th iterate 

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
   
   % compute current gradient
   gradf_xk = gradf(xk);

   % check the stopping criterion 
   if norm(gradf_xk,2) <= tol
      info = 0; % set output flag
      break;
   end
   
   % compute current hessian
   hessf_xk = hessf(xk);
   
   % try solving Newton system
   [dkN, R] = linsolve(hessf_xk, -gradf_xk);
   
   % check condition of hessian; if too bad, do not use the Newton direction
   if (R < 1.0d-8)
      dkN = 0;    % error mark: ensures that gradient direction will be used
   end
 
   % select search direction according to direction condition (7.1)
   if (norm(dkN,2) >= 1.0d-15) && ( -gradf_xk.' * dkN  >=  min(a1, a2*norm(dkN,2)^p) * norm(dkN,2)^2 )
      dk = dkN;
      fprintf('N'); % display N for newton direction
   else
      dk = - gradf_xk;
      fprintf('G'); % display G for gradient direction
   end
   
   % determine Armijo step size in direction dk
   sigmak = armijo(f, gradf_xk, xk, dk, gamma, beta, 0);

   % perform the step
   xk = xk + sigmak * dk;
   
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


return % finito

end % end of function