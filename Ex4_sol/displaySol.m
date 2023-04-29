function displaySol(fignum, fv, X, desc)
   % visualization helper
   figure(fignum); clf; hold on                 % prepare figure window
   [x1, x2]  = meshgrid(-3:0.05:3, -2:0.05:4);  % get the evaluation grid
   fvals     = fv(x1, x2);                      % evaluate the vectorized function
   hsurf     = surf(x1, x2, fvals);             % surface plot
   contours  = [0.1 1 5 10  100   1000];        % contour lines at specified f values
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
   plot3(X(1,:), X(2,:), zz, 'r.-', 'LineWidth', 2.0, 'MarkerSize', 12.0);  % plot on the z=0-plane
   plot3(X(1,:), X(2,:), fX, 'k.-', 'LineWidth', 0.5, 'MarkerSize',  8.0);  % plot on the function
      
end


