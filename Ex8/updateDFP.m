function Hnew = updateDFP(H, s, y)
   % function Hnew = updateDFP(H, s, y)
   %
   % Explicitly forms the DFP update of Hessian approximation H
   %
   % INPUT:  H --> current Hessian approximation
   %         s --> step in iterates
   %         y --> step in gradients
   %
   % OUTPUT: Hnew --> DFP update of H
   %
   
   % some intermediates for re-use
   y_Hs = y - H * s;
   yTs  = y.' * s;
   
   % DFP update
   Hnew = H  +  (y_Hs * y.' + y * y_Hs.')/(yTs) ...
             -  (y_Hs.' * s)/(yTs)^2 * y*y.';
   
end