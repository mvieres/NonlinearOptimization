function X = L_BFGS(f,gradf,x0,m,tol,maxit)
%
%
%
%
%
%
%
%
%
%
%
chunksize = 100;
x0 = reshape(x0,[],1);

X = zeros(length(x0),chunksize);
s = zeros(length(x0),chunksize);
y = zeros(length(x0),chunksize);
rho = zeros(1,100);


X(:,1) = x0;
gradf_x0 = gradf(x0);

% Stopping for x0
if norm(gradf_x0) <= tol
    return
end

for i = 1:maxit

    % Check for enough storage within arrays
    if i >= chunksize
        
        X(:,(end+1):(end+100)) = zeros(length(x0),100);
        s(:,(end+1):(end+100)) = zeros(length(x0),100);
        y(:,(end+1):(end+100)) = zeros(length(x0),100);
        rho((end+1):(end+100)) = zeros(1,100);
        chunksize = chunksize +100;
    end

    xk = X(:,i);

    % Determine search directio
    if i == 1
        d = -gradf_x0;
        gradfxk = gradf_x0;
    else
        % Chose memory size m (for k<m it has to be reduced to avoid
        % negative indicies)
        if  i<=m
            m_chosen = i-1;
        else
            m_chosen = m;
        end
        gradfxk = gradf(xk);
        % Get product of  Bk and gradfxk
        bkgrad = lbfgsinverse_evalBkd(gradfxk, i, m_chosen, rho, s, y);
        d = -bkgrad;  
    end
    
    % Stepsize strategy
    sigma = powellwolfe(f, gradf, xk, d);

    % Update the next iterate
    step = sigma*d;
    xk_plus1 = xk + step;
    X(:,i+1) = xk_plus1;

    % Breaking criterium
    gradfxk_plus1 = gradf(xk_plus1);
    if norm(gradfxk_plus1) <= tol
        break
    end

    % Compute / Store Curvature Information
    s(:,i) = step;
    y(:,i) = gradfxk_plus1 - gradfxk;
    rho(i) = 1/dot(y(:,i),s(:,i));
    
end
    X = X(:,1:i);
end