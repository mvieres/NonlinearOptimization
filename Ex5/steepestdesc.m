function [X, info] = steepestdesc(f, gradf, x0, tol, maxit, varargin)
% X = steepestdesc(f, gradf, x0, tol, maxit)
% X = steepestdesc(f, gradf, x0, tol, maxit, stepsizefunc, stepsizeparam)
% [X, info] = steepestdesc(...)
%
% Locates a stationary point of function f using the method of steepest descent.
%
% INPUT:        f --> function handle to objective (accepting one input variable)
%           gradf --> function handle to objective's gradient (accepting one input variable)
%              x0 --> initial guess (starting point)
%             tol --> stopping condition tolerance:   ||gradf|| <= tol
%           maxit --> maximum number of iterations
%    stepsizefunc --> handle to step-size control               [default: @powellwolfe]
%                     (a function func(f, gradf, x, d, params)
%   stepsizeparam --> parameter structure for stepsizefunc      [default: empty ]
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
stepsizefunc   = @powellwolfe;
stepsizeparams = struct('gamma', 1.0d-4, 'eta', 0.9);

% process optional args
if (nargin >= 6)                    % nargin = TOTAL number of input args
   stepsizefunc   = varargin{1};    % varargin{1} = first optional input arg
   stepsizeparams = varargin{2};    % varargin{2} = second optional input arg
end

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
   
   % compute gradient at xk
   gradf_xk = gradf(xk);
   
   % steepest descent search direction: negative gradient
   dk = - gradf_xk;
   
   % compute stepsize
   sigmak = stepsizefunc(f, gradf, xk, dk, stepsizeparams);
   
   % perform the step
   xk = xk + sigmak * dk;
   
   % store it in the result matrix; grow-and-copy only if necessary
   historyIdx = historyIdx + 1;
   if historyIdx > historySize
      X = [X, zeros(dimx, historyChunk)];  %#ok<AGROW>  % Suppress warning
      historySize = historySize + historyChunk; %    as we know what we do 
   end
   X(:,historyIdx) = xk;
   
end

% shrink result array to appropriate size
X = X(:, 1:historyIdx);

% finito of function
return

end % end of function