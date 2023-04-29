function [sigma] = armijo(f,gradf,x,d,gamma,beta)

% sigma = armijo(f, gradf, x, d, gamma, beta)
%
% Computes the Armijo step length for given search direction by a ...
% backtracking approach.
%
% INPUT:     
%   f --> function handle to objective (accepting one input variable)
%   gradf --> function handle to objective's gradient
%   x --> current iterate
%   d --> search direction
%   gamma --> parameter for descent sufficiency   [default: 1.0d-4 ]
%   beta --> parameter for step size shrinking   [default: 0.5    ]
%
% OUTPUT: 
%      sigma --> Armijo step size 
%                


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
    end
end