% Rosenbrock: objective and its gradient
fv = @(x1,x2) 100 * ( x2  -  x1.^2).^2 + (1 -  x1 ).^2;                               % vectorized for fast visualization

% Test function: Powell's function
% fv    = @(x1,x2,x3,x4) (x1+10*x2).^2 + 5*(x3-x4).^2 + (x2-2*x3).^4 + 10*(x1-x4).^4;

% Test function:  LOCAL vs GLOBALIZED NEWTON
% fv = @(x1,x2) log(exp(x1) + exp(-x1)) + x2.^2;

% Auto-generate f, gradf, hessf from fv using symbolic toolbox.
% Note: The generated code is usually not efficient, and contains expressions that might lead to cancellation.
% But: It's very convenient!
xsymb = num2cell(sym('x',[1 nargin(fv)]));
fvsym = fv(xsymb{:});
f     = matlabFunction(         fvsym , 'Vars', {[xsymb{:}].'});
gradf = matlabFunction(gradient(fvsym), 'Vars', {[xsymb{:}].'});
hessf = matlabFunction( hessian(fvsym), 'Vars', {[xsymb{:}].'});

                                
% initial points
clear x0;
x0{1} = [  1.0 ;  -0.5]; 
x0{2} = [ -1.2 ;   1.0]; 
x0{3} = [ 10.0 ;  12.0];
x0{4} = [ 33.0 ; -62.0];
%clear x0; x0{1} = [ 3 ; -1 ; 0 ; 1 ]; % for Powell's function


% augment function with an evaluation counter (slow, but just for the statistics)
f     = augmentWithEvalCounter(f);
gradf = augmentWithEvalCounter(gradf);
hessf = augmentWithEvalCounter(hessf);


% optimization settings
tol   = 1.0d-9;   % stopping criterion || gradf || <= tol
maxit = 100;      % maximum number of iterations 
beta  = 0.5;      % Armijo step size shrink factor
gamma = 1.0d-4;   % Armijo sufficient descent parameter
a1    = 1d-6;     %
a2    = 1d-6;     % Newton step acceptance parameters
p     = 0.1;      %

% prepare storages
k       = 0;
x0count = length(x0);
desc    = cell(x0count,1); 
X       = cell(x0count,1); 
info    = cell(x0count,1); 


% do all simulations
for i = 1:x0count  % cycle through initial guesses x0
   k = k + 1;
   f('reset');               % reset evaluation counters
   gradf('reset');           % reset evaluation counters
   hessf('reset');           % reset evaluation counters
%  [X{k}, info{k}] = localnewton(f, gradf, hessf, x0{i}, tol, maxit);
   [X{k}, info{k}] = globalnewton(f, gradf, hessf, x0{i}, tol, maxit, beta, gamma, a1, a2, p);
   fn = f('count');          % retrieve evaluation counters
   gradfn = gradf('count');  % retrieve evaluation counters
   hessfn = hessf('count');  % retrieve evaluation counters
   itern  = size(X{k},2)-1;  % minus 1 for x0
   desc{k} = sprintf('#%2d --> converged = %s:  #it = %5d', k, mat2str(~info{k}), itern);
   fprintf('\n');
   fprintf('#%3d:  x0 = [%6.2f ; %6.2f]', k, x0{i}(1), x0{i}(2));
   fprintf('   x = [%7.4f ; %7.4f]', X{k}(1,end), X{k}(2,end));
   fprintf('   f(x) = %-8.5f', f(X{k}(:,end)));
   fprintf('   #it = %-5d   #f = %-5d   #gradf = %-5d   #hessf = %-5d', itern, fn, gradfn, hessfn);
   fprintf(' \n');
end


% display solution in individual figures
if length(X{1}(:,1))~=2, return; end  % only display 2-dimensional problems
figbase = 333;
for k = 1:x0count
   if any(any(isnan(X{k}))), continue; end   % skip NaNs (for local Newton)
   displaySol(figbase+k, fv, X{k}, desc{k})
end


