%%% Sheet 7 Main

% Parameters

maxiter = 10000;
tol = 1.0d-6;
eta = 0.9;
gamma = 1.0d-4;


alpha = 20;
n = 100;
m= 20;


rng(1)
x0 = sqrt(20)*randn([n, 1]);
f = @(x) extRosenbrock(x, alpha);
gradf = @(x) extRosenbrockGradient(x, alpha);
X = L_BFGS(f,gradf,x0,m,tol,maxiter);



