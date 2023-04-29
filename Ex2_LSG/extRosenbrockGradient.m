function gradf = extRosenbrockGradient(x, alpha)
  % gradf = extRosenbrockGradient(x, alpha)
  %
  % INPUT:  x --> position to evaluate
  %     alpha --> parameter value 
  %
  % OUTPUT: gradf --> gradient value at x
  %
  % counting calls
  persistent counter
  if (nargin==1)
     if (strcmpi(x, 'reset')), counter = 0;     return; end
     if (strcmpi(x, 'count')), gradf = counter; return; end
     error('Invalid call');
  end
  counter = counter + 1;
 
  % dimension of x
  n = length(x);
  
  % indices
  idxEVEN = 2:2:n;
  idxODD  = 1:2:n;
  
  % calculate helper
  helper = 2 * alpha * (x(idxEVEN) - x(idxODD).^2);
  
  % gradient
  gradf = zeros(n, 1);
  gradf(idxODD)  = helper .* (-2 * x(idxODD)) - 2*(1 - x(idxODD));
  gradf(idxEVEN) = helper;

end