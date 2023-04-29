function status = check64(f,gradf,x,d,sigma,gamma)
    % %CHECK65 checks if condition 6.4 is stisfied
    % Syntax: status = check64(f,gradf,x,d,sigma,gamma)
    %Input:
    %   x:      Iterate
    %   d:      search direction
    %   f:      function handle
    %   gradf:  function handle gradient f
    %   sigma:  step size (scalar)
    %   gamma:  
    %
    %Output:
    %   status: Boolean returning true if condiion 6.4 is
    %   satisfied

    % Controll
    %assert(~isvector(d),'d is not a vector')
    %assert(~isscalar(sigma),'sigma is not a scalar')

    %Main
    if (f(x+sigma*d)-f(x)) > sigma*gamma*dot(gradf(x),d)
        status = false;
    else
        status = true;
    end
end