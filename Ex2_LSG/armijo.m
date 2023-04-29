function [sigma, xnew, fnew] = armijo(f, gradf, x, d, gamma, beta, varargin)
% sigma = armijo(f, gradf, x, d, gamma, beta)
% sigma = armijo(f, gradf, x, d, gamma, beta, minstep)
% [sigma, xnew, fnew] = armijo(...)
%
% Computes the Armijo step length for given search direction by a backtracking approach.
%
% INPUT:     f --> function handle to objective (accepting one input variable)
%        gradf --> function handle to objective's gradient (accepting one input variable)
%                  or a vector containing the gradient of f at x
%            x --> current iterate
%            d --> search direction
%        gamma --> parameter for descent sufficiency
%         beta --> parameter for step size shrinking
%      minstep --> minimum step size (optional, default: 1.0d-16)
%
% OUTPUT: 
%      sigma --> Armijo step size 
%                or defined minimum step size if no sufficient descent could be found.
%       xnew --> calculated new iterate  [optional]
%       fnew --> function value at xnew  [optional]
%

% set defaults
minstep = 1.0d-16;

% process optional arguments
if (nargin>6), minstep = varargin{1}; end

% if gradf is a function, evaluate it
if ~isvector(gradf)
   gradf = gradf(x);
end

% setup line search and first step
sigma   = 1;       
xnew    = x + sigma*d;
fold    = f(x);
fnew    = f(xnew);

% precompute sufficient descent
gamma_gradf_d = gamma * gradf.' * d;

% perform line search
while (fnew - fold) > sigma * gamma_gradf_d
   sigma = beta * sigma;                  % shrink step size by beta
   xnew = x + sigma*d;                    % trial step
   fnew = f(xnew);                        % function value at trial step
   % stop if minimum step size is reached
   if (sigma <= minstep)
      fprintf('Armijo: Minimum stepsize!\n');
      break; 
   end
end

% if the program reaches this point, sigma is either a step size
% with sufficient descent or the minimum step size.



end  % finito of function