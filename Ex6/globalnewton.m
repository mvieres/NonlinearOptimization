function X = globalnewton(f, gradf, hessf, x0, e, maxit, beta, gamma, a1, a2, p)
%  X = globalnewton(f, gradf, hessf, x0, e, maxit, beta, gamma, a1, a2, p)
%
% Input: 
%   f -->   Function handle
%   gradf --> Gradient of f 
%   x0 --> Starting vector
%   e --> Error Tolerance (scalar)
%   maxit --> Maximal steps of iteration 
%   beta
%   gamma
%   a1
%   a2
%   p
%
% Output:
%   X --> Matrix containing Iterates

% Init
X = zeros(size(x0,1),maxit+1);
X(:,1) = reshape(x0,[],1);
disp('Starting value: ')
disp(x0)
for i = 1:maxit
    gradf_xk = gradf(X(:,i));

    % Stopping Criterium
    if norm(gradf_xk)<= e
        break
    end

    % Compute search direction d

    % Check for singular Matrix
    hessian = hessf(X(:,i));
    if det(hessian)<= e % det terrible slow, but our matrix is 2x2, hence
        % it should work
        d = -gradf_xk;
        fprintf('G')
    else
        % Compute / Chose Newton Direction
        d = hessian\(-gradf_xk);

        % Check Condition
        norm_dN = norm(d);
        if norm_dN >= e
            if dot(-gradf_xk,d) < min([a1,a2*norm_dN^p])*norm_dN^2
                %Chose Steepest Descent Direction
                d = -gradf_xk;
                fprintf('G')
            else
                fprintf('N')
            end
        end
    end

    sigma = armijo(f,gradf,X(:,i),d,gamma,beta);

    % Output for full / reduced step
    if sigma == 1
       fprintf('F')
    else
        fprintf('r')
    end

    % Updating the iterate
    X(:,i+1) = X(:,i) + sigma*d;
end
% Cutting of the Matrix
X = X(:,1:i);

% Print Out
fprintf('\n')
disp(['Numeber iterations: ', num2str(length(X(1,:))-1)])
fprintf('\n')

end