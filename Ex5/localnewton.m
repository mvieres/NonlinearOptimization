function X = localnewton(~,gradf, hessf, x0, e, maxit)
    % LOCALNEWTON performs the Local newton method for optimization problems
    % Syntax: X = localnewton(f, gradf, hessf, x0, e, maxit)
    %
    % Inputs:
    %   f:      Function handle -- (Function to be minimized)
    %   gradf:  Function handle -- Vector (gradient of f)
    %   hessf:  Function handle -- Matrix (Hessian of f)
    %   x0:     Starting piont (Vector)
    %   e:      error tolerance
    %   maxit:  Integer (Maximal Iteration steps)
    %
    % Outputs:
    %   X:      Vector of iterates

    % Init:
    X = zeros(size(x0,1),maxit+1);
    x0 = reshape(x0,[],1);
    X(:,1) = x0;
    
%     for i = 1:maxit
%         gradf_x = gradf(X(:,i));
%         if norm(gradf_x) > e
%             %Compute newton step
%             d = linsolve(hessf(X(:,i)),-1*gradf_x);
% 
%             %Update the iterate
%             X(:,i+1) = X(:,i) + d;
%         else
%             X = X(:,1:i); % Output
%             break
%         end
%     end

    for i = 1:maxit
        gradf_x = gradf(X(:,i));
        
        if norm(gradf_x) <= e
            X = X(:,1:i);
            break
        end
        % Compute search direction
        hessf_xk = hessf(X(:,i));
        d = linsolve(hessf_xk, -gradf_x);
        %Update the iterate
        X(:,i+1) = X(:,i) + d;
    end
end