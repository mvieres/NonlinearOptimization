%%% Sheet 7 Main

% Parameters

maxiter = 10000;
tol = 1.0d-6;
eta = 0.9;
gamma = 1.0d-4;


% alpha = 20;
% 
% m= 20;
% 
% 
% rng(300)
% x0 = sqrt(20)*randn([n, 1]);
% f = @(x) extRosenbrock(x, alpha);
% gradf = @(x) extRosenbrockGradient(x, alpha);
% X = L_BFGS(f,gradf,x0,m,tol,maxiter);

alpha = [1, 20];
n = [100, 1000, 10000];
m = [1, 5, 20];

results = zeros(5,1);

tic
for i = 1:2
    for j = 1:3
        for h = 1:3
            tic
            rng(300)
            x0 = sqrt(20)*randn([n(j), 1]);
            f = @(x) extRosenbrock(x, alpha(i));
            gradf = @(x) extRosenbrockGradient(x, alpha(i));
            X = L_BFGS(f,gradf,x0,m(h),tol,maxiter);
            time = toc;
            num_fun_eval = 0;
            num_it = length(X(1,:))-1;
            
            % Abspeichern Parameter und Ergebnisse
            results(1,end) = alpha(i);
            results(2,end) = n(j);
            results(3,end) = m(h);
            results(4,end) = num_fun_eval;
            results(5,end) = time;
            results(:,end+1) = zeros(5,1);
        end
    end
end
results = results(:,1:(end-1));