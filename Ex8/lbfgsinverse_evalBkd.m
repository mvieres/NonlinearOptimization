function d = lbfgsinverse_evalBkd(gradfxk, k, m, rho, s, y)
% SYNTAX d = lbfgsinverse_evalBkd(gradfxk, k, m, rho, s, y)
% Performes a Leftinverse BFGS two loop recursion
% Input:
%   gradfxk --> gradient of f at Iterate x_k (vector)
%   k       --> Iterate step (integer)
%   m       --> memory size (integer)
%   s       --> Matrix of iterates increments (array)
%   y       --> Matrix containing gradient incremets of iterates (array)
%   rho     --> Inverse of dot-product of s_k and y_k (1xk array) 
%
% Output:
%   d --> Product of B_k and gradfxk (array)

q = gradfxk;
alpha = zeros(1,m);
for i = (k-1):-1:(k-m)
    alpha(k-i) = rho(i)*dot(s(:,i),q);
    q = q - alpha(k-i)*y(:,i);
end

d = (dot(s(:,k-1),y(:,k-1))/dot(y(:,k-1),y(:,k-1)))*q;

for i = (k-m):(k-1)
    beta = rho(i)*dot(y(:,i),d);
    d = d + s(:,i)*(alpha(k-i)-beta);
end

end