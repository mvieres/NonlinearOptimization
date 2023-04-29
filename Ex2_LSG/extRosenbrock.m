function f = extRosenbrock(x, alpha)
  % f = extRosenbrock(x, alpha)
  %
  % INPUT:  x --> position to evaluate
  %     alpha --> parameter value 
  %
  % OUTPUT: f --> function value
  %
  
  % counting calls
  persistent counter
  if (nargin==1)
     if (strcmpi(x, 'reset')), counter = 0; return; end
     if (strcmpi(x, 'count')), f = counter; return; end
     error('Invalid call');
  end
  counter = counter + 1;
 
  % dimension of x
  n = length(x);
  
  % indices
  idxEVEN = 2:2:n;
  idxODD  = 1:2:n;
  
  % calculate
  f = sum(  alpha * ( x(idxEVEN) - x(idxODD).^2 ).^2   +   (1 - x(idxODD)).^2 );
  
end