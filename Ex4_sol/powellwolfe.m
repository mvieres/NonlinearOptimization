function sigma = powellwolfe(f, gradf, x, d, varargin)
% sigma = powellwolfe(f, gradf, x, d)
% sigma = powellwolfe(f, gradf, x, d, params)
% sigma = powellwolfe(f, gradf, x, d, gamma, eta)
%
% Determines a steplength sigma that fulfills the Wolfe conditions.
%
% INPUT:  f --> function handle to the objective
%     gradf --> function handle to the objective's gradient
%         x --> current iterate
%         d --> search direction
%     gamma --> sufficient descent parameter  [default: 1d-4]
%       eta --> initial slope parameter       [default: 0.9 ]
%    params --> structure with fields
%               .gamma --> sufficient descent parameter
%               .eta   --> initial slope parameter
%    
%
% OUTPUT:  sigma --> steplength fulfilling the Wolfe conditions.
%

% set defaults
gamma = 1.0d-4;
eta   = 0.9;

% process optional input
if (nargin == 5)      % user has specified a single parameter structure
   params = varargin{1}; 
   gamma  = params.gamma;
   eta    = params.eta;
elseif (nargin == 6)  % user has specified two inputs
   gamma = varargin{1};
   eta   = varargin{2}; 
end


% Notation:    xold = x^k        = x
%              xnew = x^{k+1}    = x^k + sigma * d^k
%              fold = f(x^k)     = f(x)
%              fnew = f(x^{k+1}) = f(x^k + sigma*d^k)

% precompute some intermediates
gradfold      = gradf(x);         % gradf(xk) = old/current gradient
fold          = f(x);             % f(xk)     = old/current function value
gradf_d       = gradfold.' * d;   % gradf(xk)'*d
gamma_gradf_d = gamma * gradf_d;  % gamma * gradf(xk)'*d
eta_gradf_d   = eta * gradf_d;    % eta * gradf(xk)'*d

% used variables:
left  = 1;  % init:    left = sigma- in the lecture
right = 1;  % init:   right = sigma+ in the lecture
            % noinit:   mid = sigma  in the lecture



% Interval Phase: 
% ===============
% Determine a "left" value that satisfies sufficient descent
% and a "right" value that does not. 

if satisfies_sufficient_descent(left)
   if satisfies_curvature_condition(left)
      sigma = left;
      return  % STOP with sigma = left
   end
   % determine smallest value 2^k (k=1,2,3,...)
   % such that the sufficient descent condition IS NOT fulfilled
   for k = 1:1023                % note: 2^1024 == inf on most 64-bit machines
      right = 2 * right;
      if ~(satisfies_sufficient_descent(right))  % ~ = NOT
         break;                  % exit the for loop
      end
   end
   left = right/2;               % adjust left interval limit
else
   % determine largest value 2^-k (k=1,2,3,...)
   % such that the sufficient descent condition IS fulfilled
   for k = 1:1074                % note: 2^-1075 == 0 on most 64-bit machines
      left = left / 2;
      if satisfies_sufficient_descent(left)
         break;                  % exit the for loop
      end
   end
   right = left*2;               % adjust right interval limit
end

% at this point, "left" does satisfy the sufficient descent condition
% and "right" does not satisfy it.

% Bisection Phase: 
% Shrink the interval [left, right] until a Wolfe step size is determined
while ~(satisfies_curvature_condition(left))    % ~ = NOT
   mid  = (left + right) * 0.5;
   if satisfies_sufficient_descent(mid)
      left  = mid;   % shrink: move left limit to mid
   else
      right = mid;   % shrink: move right limit to mid
   end
end

% finito: "left" fulfills both Wolfe conditions
sigma = left;

return




   % Helpers:  These *nested functions* have access to the variables 
   % ========  lying in the workspace of the powellwolfe function
   %           Note: This is primary for readability of the code.
   %
   
   function flag = satisfies_sufficient_descent(sigma)
      % (6.2): sufficient descent condition
      % f(xnew)  <=  f(x) + sigma * gamma * gradf(x).' * d
      xnew = x + sigma * d;  % NOTE: inefficient recomputation!
      fnew = f(xnew);
      flag = ( fnew <= fold + sigma * gamma_gradf_d );
   end

   function flag = satisfies_curvature_condition(sigma)
      % (6.3): curvature condition
      % gradf(xnew).'d  >=  eta * gradf(x).' * d
      xnew = x + sigma * d;  % NOTE: inefficient recomputation!
      gradfnew = gradf(xnew); 
      flag = ( gradfnew.' * d >= eta_gradf_d );
   end


end