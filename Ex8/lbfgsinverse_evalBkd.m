function d = lbfgsinverse_evalBkd(gradfxk, k, m, rho, s, y)
% SYNTAX d = lbfgsinverse_evalBkd(gradfxk, k, m, rho, s, y)
% Performes a Leftinverse BFGS two loop recursion
% Input:
%   gradfxk --> gradient of f at Iterate x_k (vector)
%   k --> Iterate step (integer)
%   m --> memory size (integer)
%   rho --> 
%   s -->
%   y -->
%
% Output:
%   d --> Product of B_k and gradfxk

q = gradfxk;
alpha = zeros(1,m);
for i = (k-1):-1:max(k-m,1)
    alpha(k-i) = rho(i)*dot(s(:,i),q);
    q = q - alpha(k-i)*y(:,i);
end

d = (dot(s(:,k-1),y(:,k-1))/dot(y(:,k-1),y(:,k-1)))*q;

for i = (k-m):(k-1)
    beta = rho(i)*dot(y(:,i),d);
    d = d + s(:,i)*(alpha(k-i)-beta);
end

end