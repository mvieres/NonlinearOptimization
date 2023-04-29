function status = check65(gradf,x,d,sigma,eta)
    %CHECK65 checks if condition 6.5 is stisfied
    % status = check65(f,gradf,x,d,sigma,eta)
    %Input:
    %   x:      Iterate
    %   d:      search direction
    %   gradf:  function handle gradient f
    %   sigma:  step size (scalar) 
    %   eta:    
    %
    %Output:
    %   status: Boolean returning true if condiion 6.5 is
    %   satisfied, flase otherwise

    % Controll
    %assert(~isvector(d),'d is not a vector')
    %assert(~isscalar(sigma),'sigma is not a scalar')
    
    if dot(gradf(x+sigma*d),d) < eta*dot(gradf(x),d)
        status = false;
    else
        status = true;
    end
end