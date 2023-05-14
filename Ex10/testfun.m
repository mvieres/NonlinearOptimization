diary("output_console.txt")
diary("on")
% setup functions
f = cell(4,1);                                                             % functions:
f{1} = @(x1,x2,x3,x4,x5) x1.^2 + x2.^2 + x3.^2 + x4.^2 + x5.^2 - 10;       % 1: constraint h1
f{2} = @(x1,x2,x3,x4,x5) x2.*x3 - 5*x4.*x5;                                % 2: constraint h2
f{3} = @(x1,x2,x3,x4,x5) x1.^3 + x2.^3 + 1;                                % 3: constraint h3
f{4} = @(x1,x2,x3,x4,x5) exp(x1*x2*x3*x4*x5) - 0.5*(x1.^3 + x2^3 +1).^2;   % 4: objective
n = length(f);

% auto-generate f, gradf, hessf using symbolic notations
fprintf('Creating derivatives..')
funs = cell(n,1); grads = cell(n,1); hesss = cell(n,1);
symfvs = cell(n,1); symgrads = cell(n,1); symhesss = cell(n,1);
for i = 1:length(f)
   symx  = num2cell(sym('x',[1 nargin(f{i})]));
   symfvs{i}   = f{i}(symx{:});
   symgrads{i} = gradient(symfvs{i},[symx{:}]);
   symhesss{i} = hessian(symfvs{i},[symx{:}]);
   funs{i}  = matlabFunction(symfvs{i}  , 'Vars', {[symx{:}].'});
   grads{i} = matlabFunction(symgrads{i}, 'Vars', {[symx{:}].'});
   hesss{i} = matlabFunction(symhesss{i}, 'Vars', {[symx{:}].'});
   fprintf('.');
end
fprintf('Done!\n');

% build Lagrangian gradient and hessian
% index 1-3:  h (constraints)   last index 4: f (objective)
gradL = @(x,mu)  grads{4}(x)  +  mu(1)*grads{1}(x) + mu(2)*grads{2}(x) + mu(3)*grads{3}(x);
hessL = @(x,mu)  hesss{4}(x)  +  mu(1)*hesss{1}(x) + mu(2)*hesss{2}(x) + mu(3)*hesss{3}(x);
% for display and comparison: build symbolic versions
symmu = num2cell(sym('mu',[1 n-1]));
symgradL = symgrads{4}  +  symmu{1}*symgrads{1} + symmu{2}*symgrads{2} + symmu{3}*symgrads{3};
symhessL = symhesss{4}  +  symmu{1}*symhesss{1} + symmu{2}*symhesss{2} + symmu{3}*symhesss{3};

% assemble constraint function h and Jacobian of constraints Jh
h  = @(x) [funs{1}(x)    ; funs{2}(x)    ; funs{3}(x)   ];
Jh = @(x) [grads{1}(x).' ; grads{2}(x).' ; grads{3}(x).'];

% as a remark: function creation like this (by symbolic computations) is inefficient!
% manual implementation of the respective functions would tremendously increase performance
% however, the implementational effort would be much higher
% also: the nesting of anonymous functions decreases performance

% initial guesses
xs  = [-1.7171, 1.5957, 1.8272, -0.76364, -0.76364].'; % solution
x0  = [-1.80, 1.40, 1.90, -0.80, -0.80].';  
mu0 = zeros(3,1);

% tolerance and maxit
tol   = 1.0d-9;
maxit = 100;

% call to optimizer
disp('Starting localSQP...')
tic
[X, MU, info] = localSQP(gradL, hessL, h, Jh, x0, mu0, tol, maxit);
toc
disp('localSQP finished.')
disp('')

% comparison to fmincon
ff = funs{4};
hh = @(x) deal([], h(x));
opts = optimoptions('fmincon','Algorithm','sqp','ConstraintTolerance',tol,'OptimalityTolerance',tol);
disp('Starting fmincon...')
tic
[xopt,fval,exitflag,output,lambda,gradvec,hessmat] = fmincon(ff, x0, [], [], [], [], [], [], hh, opts);
toc
disp('fmincon finished.')
disp('')

% solution comparison
fprintf('Solution Comparison:\n  localSQP   fmincon   solution\n')
disp([X(:,end), xopt(:), xs(:)]);
fprintf('Multiplier Comparison:\n  localSQP   fmincon\n')
disp([MU(:,end), lambda.eqnonlin]);

diary("off")