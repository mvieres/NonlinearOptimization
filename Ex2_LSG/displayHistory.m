function displayHistory(history, fieldwidth, acc, kmin, kmax)

% defaults
if (nargin < 2), fieldwidth = 10;                end
if (nargin < 3),        acc =  6;                end
if (nargin < 4),       kmin =  1;                end
if (nargin < 5),       kmax = length(history)-1; end

% extract the data from history
xk     = [history.xk]; 
fk     = [history.fk]; 
dk     = [history.dk]; 
sigmak = [history.sigmak];

% generate generic format strings
fmt_int    = sprintf('| %%%dd ', fieldwidth);
fmt_dble   = sprintf('| %%%d.%df ', fieldwidth, acc);
fmt_header = sprintf('| %%%ds ', fieldwidth);

% determine the number of components of xk and prepare the format for it
dim_xk = min(5, length(xk(:,1))); % display no more than the first 5 components

% display the header once
fprintf(fmt_header, 'k');
for l = 1:dim_xk
   fprintf(fmt_header, sprintf('xk_%d', l));
end
fprintf(fmt_header, 'fk');
fprintf(fmt_header, 'sigmak');
fprintf(fmt_header, '||dk||');
fprintf('\n');


% display the iterates
for k = max(1,kmin) : min(length(history),kmax)
   displayIterate(k);
end

% display last iterate
if kmax < length(history)
   fprintf('%10s\n',  '. . .');
   displayIterate(length(history));
end

% finito
return



% helper for displaying the k-th iterate
function displayIterate(k)
   fprintf(fmt_int, k);                 % iteration counter
   for i = 1:dim_xk                     
      fprintf(fmt_dble, xk(i,k));       % components of iterate xk
   end
   fprintf(fmt_dble, fk(k));            % fk
   fprintf(fmt_dble, sigmak(k));        % sigmak
   fprintf(fmt_dble, norm(dk(:,k),2));  % 2-norm of dk
   fprintf('\n');                       % emit new line   
end



end


