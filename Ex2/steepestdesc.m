function X = steepestdesc(f,gradf,x0,e,maxit)
    %Init
    x_k = x0; % should be a column vector
    i = 1; %counting each Iteration step starting with 1 sinc matlab starts with 1
    X = x0;

    %Parameter for armijo
    gamma = 10^-3;
    beta = 0.5;
    %Loop
    while norm(gradf(x_k)) > e
        if i <= maxit
            d = -1*gradf(x_k); %Computing the gradient
            sigma = armijo(f,gradf,x_k,d,gamma,beta); %Sigma from armijo
            x_k = x_k + sigma*d; %Computing the new iteration value
            X(:,end+1) = x_k; %Saving the new iterate in each column
            i = i + 1; %updating the iteration step
            %if (length(X(:,1)) < 2) && (i < min(6, length(X(1,:))))
                diary Ex2_Text
                disp(['Step: ', num2str(i-1)])
                disp(['Step Size: ', num2str(sigma),])
                disp(['Iterative:', num2str(x_k)])
                disp(['Search direction: ', num2str(d)])
                disp(['Function value: ', num2str(f(x_k))])
                diary off
            %end
        else
           break
        end
    end
end