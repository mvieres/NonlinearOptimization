function sigma = powellwolfe(f, gradf, x, d, gamma, eta)
    % POWELLWOLFE performes a line search based on Powell and Wolfe
    % Syntax: sigma = powellwolfe(f, gradf, x, d, gamma, eta)
    %
    % Inputs:
    %   f:       function handle
    %   gradf:   function handle (Gradient of f)
    %   d:       Search direction (n-dimensional vector)
    %   gamma:   float in (0,0.5)
    %   eta:     float in (gamma,1)
    %
    % Output: 
    %   sigma:   stepsize
    
    % Init
        sigma_minus = 1.5;
        N = 1000; % Reducing the set with infinite many elements from Algorithmus

    % Controll
    assert(isvector(x),'x is not a vector')
    assert(isvector(d),'d is not a vector')

    % If loop for: sigma minus satisfies Armijo codition 6.4
    status = check64(f,gradf,x,d,sigma_minus,gamma);
    if status
        %If loop for sigma_minus satisfying also Curvature condition 6.5 
        status = check65(gradf,x,d,sigma_minus,eta);
        if status
            %Output:
            sigma = sigma_minus;
            return
        end

        %Find smallest sigma_plus s.t. 6.4 is not satisfied
        
        sigma_plus = 2.^(1:N);
        status = arrayfun(@(sigma_plus) check64(f,gradf,x,d,sigma_plus,gamma),sigma_plus);
        values = find(status == false);
        assert(~isempty(values),'Cond 6.4 is not satisfied for all i=1,...,N')
        sigma_minus = (2^values(1))/2;
    else
        % Find smallest sigma_minus st 6.4 holds
        sigma_minus = 2.^-(1:N);
        status = arrayfun(@(sigma_minus) check64(f,gradf,x,d,sigma_minus,gamma),sigma_minus);
        values = find(status == true);
        assert(~isempty(values),'Cond 6.4 is not satisfied for all i=1,...,N')
        sigma_minus = 2^(-values(1));
        sigma_plus = 2*sigma_minus;
    end

    % While loop as long as sigma_minus does not satisfy 6.5
    status = check65(gradf,x,d,sigma_minus,eta);
    while ~status %sigma_minus does not satisfy 6.5
        sigma = (sigma_plus + sigma_minus)/2;
        status = check64(f,gradf,x,d,sigma,gamma);
        %If loop if sigma satisfies 6.4
        if status
            sigma_minus = sigma;
        else
            sigma_plus = sigma;
        end
        %Checking if updated sigma_minus satisfies 6.5
        status = check65(gradf,x,d,sigma_minus,eta);
    end
    % Output
    sigma = sigma_minus;
end