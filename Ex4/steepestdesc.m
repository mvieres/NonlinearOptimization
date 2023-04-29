function X = steepestdesc(f,gradf,x0,e,maxit)
    % STEEPESTDESC performs the steepest decent methods to compute local
    % minimum
    % Syntax: X = steepestdesc(f,gradf,x0,e,maxit)
    %
    % Input:
    %       f       --> Function handle of function f
    %       gradf   --> Function handle of gradient of f
    %       x0      --> Starting point of itereation
    %       e       --> Tolerance
    %       maxit   --> maxium of iteration steps
    % Output:
    %       X       --> Vector of Iterates

    %Init
    x_k = x0; 
    i = 2; 
    X = zeros(2,maxit);
    X(:,1) = x0;

    %Parameter for powell wolfe
    gamma = 1e-4;
    eta = 0.9;

    %
    while i <= (maxit-1)
        gradf_xk = gradf(x_k);
        if norm(gradf_xk,2) > e
            d = -gradf_xk;
            sigma = powellwolfe(f, gradf,  x_k, d, gamma, eta);
            x_k = x_k + sigma*d; 
            X(:,i) = x_k;  
        else
            break
        end
       i = i + 1;
    end
    X = X(:,1:(i-1)); 
end