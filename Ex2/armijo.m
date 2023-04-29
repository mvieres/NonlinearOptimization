function sigma = armijo(f,gradf,x,d,gamma,beta)
    %Init
    sigma = 1;
    %Armijo
    while (f(x + sigma*d) - f(x)) > (sigma*gamma*gradf(x)'*d)    
       sigma = beta*sigma; %Updating sigma
    end
end