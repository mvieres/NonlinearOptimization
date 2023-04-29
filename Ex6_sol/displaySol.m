function displaySol(fignum, fv, X, desc)

   % determine axes' bounds and make grid
   x1min = min(X(1,:));  x1max = max(X(1,:));
   x2min = min(X(2,:));  x2max = max(X(2,:));
   nsteps = 250;
   x1steps = (x1max-x1min)/nsteps;
   x2steps = (x2max-x2min)/nsteps;
   x1grid = (x1min-1):x1steps:(x1max+1);
   x2grid = (x2min-1):x2steps:(x2max+1);
   
   % evaluate and determine contour levels
   [x1, x2]  = meshgrid(x1grid, x2grid);        % evaluation grid at interesting area 
   fvals     = fv(x1, x2);                      % evaluate the vectorized function
   fmin      = min(min(fvals));                 % minimum f value
   fmax      = max(max(fvals));                 % maximum f value
   contours  = logspace(log10(fmin), log10(fmax), 8);  % contour lines at log-spaced intervals
   
   % visualization
   figure(fignum); clf; hold on                 % prepare figure window
   hsurf     = surf(x1, x2, fvals);             % surface plot
   [~,hcont] = contour(x1, x2, fvals, contours);% contour plot
   hsurf.FaceAlpha = 0.5;                       % transparency of surface plot faces
   hsurf.EdgeAlpha = 0.1;                       % transparency of surface plot edges
   hcont.LineWidth = 1.5;                       % line width of contour plot
   xlabel('x1');                                % set label for x-axis
   ylabel('x2');                                % set label for y-axis
   title(desc, 'Fontsize', 8);                  % set title
   
   % add the optimization paths
   fX = fv(X(1,:), X(2,:)); % evaluate the function at specified points
   zz = zeros(size(fX));    % for plotting in the z=0 plane
   plot3(X(1,:), X(2,:), fX, 'k.-', 'LineWidth', 1.5, 'MarkerSize', 8.0);  % plot on the function
   plot3(X(1,:), X(2,:), zz, 'r.-', 'LineWidth', 1.5, 'MarkerSize', 8.0);  % plot on the z=0-plane
   
end


