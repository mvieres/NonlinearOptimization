function X = globalnewtonlike(f, gradf, Mkf, x0, e, maxit, beta, gamma, a1, a2, p)
% Syntax X = globalnewtonlike(f, gradf, Mkf, x0, e, maxit, beta, gamma, a1, a2, p)
% GLOBALNEWTONLIKE performs the global newton like descent method 
%
% Input:
%   f --> function f (function handle)
%   gradf --> gradient of f (function handle)
%   Mkf --> Matrix f (function handle OR matrix)
%   x0 --> Starting iterate (vector)
%   e --> Null-Tolerance (Number)
%   maxit --> Maximal iterates (integer)
%   beta --> Beta parameter (compare with Skript)
%   gamma --> gamma parameter
%   a1 --> a1 parameter
%   a2 --> a2 parameter
%   p --> p parameter
%
% Output:
%   X --> Matrix with all iterates
%   Textoutput in Console


% Init
X = zeros(size(x0,1),100);
X(:,1) = reshape(x0,[],1);
chunksize = 100;
stepsizes = zeros(1,100);
step_chosen = strings([1,100]);
%Checking for constant input in Mkf
flag = isa(Mkf,'function_handle');

if ~flag
    Mkf_xk = Mkf; % constant Matrix Mkf for every iteration
end

for i = 1:maxit
    gradf_xk = gradf(X(:,i));
    neg_gradf = -gradf_xk;
    norm_gradf = norm(gradf_xk);

    %Stopping Criterium
    if norm_gradf <= e
        break
    end
    
    
    if flag
        Mkf_xk = Mkf(X(:,i));
    end

    %
    [d, R] = linsolve(Mkf_xk,neg_gradf);
    step_chosen(i) = 'N';
    if R < e %Checking for Singularity of Mkf
        d = neg_gradf;
        step_chosen(i) = 'G';
    else
        % Check Condition
        norm_dN = norm(d);
        if norm_dN >= e
            if dot(neg_gradf,d) < min([a1,a2*(norm_dN^p)])*norm_dN^2
                %Chose Steepest Descent Direction
                d = neg_gradf;
                step_chosen(i) = 'G';
            end
        end
    end

    % Compute Step Size via Armijo Step Size Strategy
    sigma = armijo(f,gradf,X(:,i),d,gamma,beta);
    stepsizes(i) = sigma;
    % Check chunksize / extend matrix X every 100 iterations
    if i >= chunksize
        chunksize = chunksize + 100;
        X = [X, zeros(size(x0,1),100)];
        stepsizes = [stepsizes, zeros(1,100)];
        step_chosen(end+1:end+100) = strings([1,100]);
    end


    % Updating the iterate
    X(:,i+1) = X(:,i) + sigma*d;

    
end
% Cutting of the Matrix
X = X(:,1:i);


% Print Out in each iteration step
        i = 1;
        norm_gradf = norm(gradf(X(:,i)));
        fprintf(' Number iteration step = %11.5g    ', i)
        fprintf('   f(x) = %11.5g    ', f(X(:,i)))
        fprintf('   normgradf = %11.5g    ', norm_gradf)
        fprintf('   xk = [%11.5g, %11.5g]   ', X(1,i), X(2,i))
        fprintf('   stepsize = %11.5g   ', stepsizes(i))    
        fprintf('   Step = %s   ', step_chosen(i))
        fprintf(' \n')

    for i = 2:(length(X(1,:))-1)  
        contraction = (norm(X(:,i+1)-X(:,i)))/(norm(X(:,i)-X(:,i-1)));
        norm_gradf = norm(gradf(X(:,i)));
        fprintf(' Number iteration step = %11.5g    ', i)
        fprintf('   f(x) = %11.5g    ', f(X(:,i)))
        fprintf('   normgradf = %11.5g    ', norm_gradf)
        fprintf('   xk = [%11.5g, %11.5g]   ', X(1,i), X(2,i))
        fprintf('   stepsize = %11.5g   ', stepsizes(i))
        fprintf('   contractionrate = %11.5g   ', contraction)
        fprintf('   Step = %s    ', step_chosen(i))
        fprintf('\n')
    end

        i = length(X(1,:));
        norm_gradf = norm(gradf(X(:,i)));
        fprintf(' Number iteration step = %11.5g   ', i)
        fprintf('   f(x) = %11.5g    ', f(X(:,i)))
        fprintf('   normgradf = %11.5g   ', norm_gradf)
        fprintf('   xk = [%11.5g, %11.5g]   ', X(1,i), X(2,i))   

        if flag
            eigen = eig(Mkf(X(:,end)));
            fprintf('   Eigenvalues = [%11.5g, %11.5g]', eigen(1), eigen(2))
        end
        fprintf(' \n')

end



