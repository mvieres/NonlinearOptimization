function sigma = armijo(f, gradf, x, d, varargin)
% sigma = armijo(f, gradf, x, d, params)
% sigma = armijo(f, gradf, x, d, gamma, beta)
%
% Computes the Armijo step length for given search direction by a backtracking approach.
%
% INPUT:     f --> function handle to objective (accepting one input variable)
%        gradf --> function handle to objective's gradient (accepting one input variable)
%                  or a vector containing the gradient of f at x
%            x --> current iterate
%            d --> search direction
%        gamma --> parameter for descent sufficiency   [default: 1.0d-4 ]
%         beta --> parameter for step size shrinking   [default: 0.5    ]
%       params --> parameter structure with fields:
%                  .gamma    as above
%                  .beta     as above
%                  .minstep  minimum step size         [default: 1.0d-100]
%
% OUTPUT: 
%      sigma --> Armijo step size 
%                or defined minimum step size if no sufficient descent could be found.
%

% set defaults
gamma = 1.0d-4;
beta  = 0.5;
minstep = 1.0d-100;  % ensure exist

% process optional input
if (nargin == 5)      % user has specified a single parameter structure
   params  = varargin{1}; 
   gamma   = params.gamma;
   beta    = params.beta;
   minstep = params.minstep;
elseif (nargin == 6)  % user has specified two inputs
   gamma = varargin{1};
   beta  = varargin{2}; 
end

% setup line search and first step
sigma   = 1;       
xnew    = x + sigma*d;
fold    = f(x);
fnew    = f(xnew);

% precompute some quantities
gradfx        = gradf(x);
gamma_gradf_d = gamma * gradfx.' * d;

% perform line search
while (fnew - fold) > sigma * gamma_gradf_d
   sigma = beta * sigma;                  % shrink step size by beta
   xnew = x + sigma*d;                    % trial step
   fnew = f(xnew);                        % function value at trial step
   % stop if minimum step size is reached
   if (sigma <= minstep)
      fprintf('Armijo: Minimum stepsize exceeded!\n');
      break; 
   end
end

% if the program reaches this point, sigma is either a step size
% with sufficient descent or the minimum step size.



end  % finito of function