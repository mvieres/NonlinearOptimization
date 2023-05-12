%%% Sheet 8 Main
diary('output_console.txt')
% Parameters

maxiter = 10000;
tol = 1.0d-6;
eta = 0.9;
gamma = 1.0d-4;

alpha = [1, 20];
n = [100, 1000, 10000];
m = [1, 5, 20];

results = zeros(5,1);

for i = 1:2
    for j = 1:3
        for h = 1:3

            f     = augmentWithEvalCounter(@(x) extRosenbrock(x, alpha(i)));
            gradf = augmentWithEvalCounter(@(x) extRosenbrockGradient(x, alpha(i)));
            
            % Reset counter
            f('reset'); 
            gradf('reset');

            % Start of calculations
            
            rng(300)
            x0 = sqrt(20)*randn([n(j), 1]);
            tic;
            %f = @(x) extRosenbrock(x, alpha(i));
            %gradf = @(x) extRosenbrockGradient(x, alpha(i));
            X = L_BFGS(f,gradf,x0,m(h),tol,maxiter);
            time = toc;
            num_it = length(X(1,:))-1;
            num_f_eval = f('count'); 
            num_gradf_eval = gradf('count');

            % Abspeichern Parameter und Ergebnisse

            results(1,end) = alpha(i);
            results(2,end) = n(j);
            results(3,end) = m(h);
            results(4,end) = num_f_eval;
            results(5,end) = time;
            results(:,end+1) = zeros(5,1);


            % Printout
            fprintf('Parameters:   dim=%5d   -   alpha=%2d   -   m=%2d', n(j), alpha(i), m(h));
            fprintf('\n')

            fprintf(['Optimization ended: iterations=%d   -   evaluations f=%d   -   evaluations gradf=%d  ...' ...
                '  -   runtime=%d'], num_it, num_f_eval, num_gradf_eval, time)
            fprintf('\n \n')

        end
    end
end
results = results(:,1:(end-1));
diary("off")