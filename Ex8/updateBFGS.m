function Hnew = updateBFGS(H, s, y)
   % function Hnew = updateBFGS(H, s, y)
   %
   % Explicitly forms the BFGS update of Hessian approximation H
   %
   % INPUT:  H --> current Hessian approximation
   %         s --> step in iterates
   %         y --> step in gradients
   %
   % OUTPUT: Hnew --> DFP update of H
   %
   
   % some intermediates for re-use
   Hs = H * s;
   yTs  = y.' * s;
   
   % BFGS update
   Hnew = H  +  (y*y.') / yTs  -  (Hs*Hs.') / (s.'*Hs);
             
   
end