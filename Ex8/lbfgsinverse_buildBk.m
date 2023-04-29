function Bk = lbfgsinverse_buildBk(k, m, B0, s, y)
   % Bk = lbfgsinverse_buildBk(k, m, B0, s, y)
   %
   % Explicity build of Bk from QN history.
   %
   % INPUT: k --> current iterate number
   %        m --> memory size
   %       B0 --> initial matrix 
   %              If B0 is empty, a unit matrix is used
   %              If B0 is scalar, a scaled unit matrix is used
   %        s --> step in iterates (cell array)
   %        y --> step in gradients (cell array)
   %
   % OUTPUT: Bk --> inverse Hessian approximation
   %
   
   % get problem size
   n = length(s{1});

   % check if B0 is specified, otherwise set it to defaults
   if isempty(B0),  B0 = eye(n);               end
   if isscalar(B0), B0 = diag(B0 * ones(n,1)); end 
      
   % restrict memory length
   m = min(m, k);
   
   % build inverse BFGS matrix by successive DFP updates
   Bk = B0;
   for i = (k-m):k-1
      Bk = updateDFP(Bk, y{i}, s{i});
   end
   
end