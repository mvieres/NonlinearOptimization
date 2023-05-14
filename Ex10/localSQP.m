function [X, MU, info] = localSQP(gradL, hessL, h, Jh, x0, mu0, tol, maxit)
% SYNTAX: [X, MU, info] = localSQP(gradL, hessL, h, Jh, x0, mu0, tol, maxit
% Inputs: 
%   gradL   --> gradient of L for gien x and mu
%   hessL   --> Hessian for L at given x and mu 
%   h       --> Function value h
%   Jh      --> Jacobian of j
%   x0      --> Initial guess
%   mu0     --> Initial guess
%   tol     --> Tolerance  
%   maxit   --> Total number iterations
%
% Outputs: 
%   X       --> Matrix
%   MU      --> Matrix containing all Lagrangian multipliers mu_k as columns
%   info    --> flag 0 = KKT point, -1 otherwise
%

% Init
x0 = reshape(x0, [],1);
mu0 = reshape(mu0, [],1);

n = length(x0);
p = length(mu0);

% Preallocation
chunksize = 100;
X = zeros(n,100);
MU = zeros(p,100);
info = -1;

X(:,1) = x0;
MU(:,1) = mu0;


for i = 1:maxit
    xk = X(:,i);
    mu_k = MU(:,i);

    h_xk = h(xk);
    gradL_k = gradL(xk,mu_k);
    grad_h = Jh(xk);
    HessL_k = hessL(xk,mu_k);

    % Stopping criterium
    if norm(h_xk) <= tol
        if norm(gradL_k) <= tol
            info = 0;
            break
        end
    end

    % Compute d_k 
    M = [HessL_k, grad_h' ; grad_h , zeros(p,p)]; % Switching transpose cause grad_h is 1xn vector
    vec = [-gradL_k ;  -h_xk];
    
    d_joint_k = M\vec;

    d_xk = d_joint_k(1:n);
    d_muk = d_joint_k(n+1:end);

    % Chunksize extension
    if i >= chunksize 
        X(:,(end+1):(end+100)) = zeros(n,100);
        MU(:,(end+1):(end+100)) = zeros(p,100);
        chunksize = chunksize +100;
    end
    
    % Updating Iterates
    X(:,i+1) = xk + d_xk;
    MU(:,i+1) = mu_k + d_muk;
end

X = X(:,1:i);
MU = MU(:,1:i);
end
