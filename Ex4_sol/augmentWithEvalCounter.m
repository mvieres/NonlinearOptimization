function fh = augmentWithEvalCounter(funhandle)
% fh = augmentWithEvalCounter(funhandle)
%
% Adds an evaluation counter to specified function handle.
%
% INPUT:   funhandle --> function to augment with evaluation counter
%
% OUTPUT:         fh --> augmented function
%
% NOTE:  This implementation is limited to function handles with
%        exactly 1 input argument and 1 output argument.
%

numEvals = 0;
fh = @f_with_eval_counter;
return;


   function y = f_with_eval_counter(x)
      
      % augment the function
      if (nargin == 1) && ischar(x)
         switch lower(x)
            case 'reset' ,  numEvals = 0;
            case 'count' ,  y = numEvals;
            otherwise    ,  error('Don''t understand: %s', x);
         end
         return
      end

      numEvals = numEvals + 1;   % increase counter
      y = funhandle(x);          % call the original function
      
   end
   
  
   
end


