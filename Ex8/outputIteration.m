function outputIteration(k, X, sigmak, gradfxk, f)
   % outputIteration(k, X, sigmak, gradfxk, f) 
   %
   % Displays information during optimization.
   %
   % INPUT:  k --> iteration number
   %         f --> function handle to objective f
   %         X --> history storage
   %    sigmak --> step size
   %   gradfxk --> gradient at xk
   %
   % OUTPUT: (none)
   %
   
   % compute/set some intermediate
   xk = X(:, k);                                     % xk   = current iterate
   if k>1, xkm1 = X(:, k-1); else, xkm1 = NaN; end   % xkm1 = previous iterate
   if k>2, xkm2 = X(:, k-2); else, xkm2 = NaN; end   % xkm2 = penultimate iterate
   progress = norm(xk-xkm1,2) / norm(xkm1-xkm2,2);
   n = length(xk);
   
   % output
   fprintf('#%3d: ', k);                           % iteration
   fprintf('| ||gradf||=%11.5g ', norm(gradfxk));  % norm of gradient
   fprintf('| crate=%11.5g ', progress);           % progress measure: contraction rate
   fprintf('| f(x)=%11.5g ', f(xk));               % function value
   fprintf('| step=%11.5g   ', sigmak);            % step size
   for i = 1:min(10,n)
      fprintf('| x(%03d)=%11.5g ', i, xk(i));      % current iterate (10 entries at max)
   end
   fprintf('\n');  % new line
end
