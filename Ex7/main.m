%%% Exercise sheet 7
diary('output_console.txt')

% Parameters
epsilon = 1e-6;
maxiter = 100;
beta = 0.5;
gamma = 1e-04;
a1 = 1e-06;
a2 = 1e-06;
p = 0.1;

% Starting iterates
x_a = [1; -0.5];
x_b = [-1.2; 1.0];
x_c = [5; 5];

%% Function 1
f1 = @(x) log(exp(x(1)) + exp(-x(1))) + x(2)^2;
fv1 = @(x1,x2) log(exp(x1) + exp(-x2)) + x2.^2;
gradf1 = @(x)[ (exp(x(1))-exp(-x(1)))*(exp(x(1)) + exp(-x(1)))^-1; 2*x(2)];
hessf1 = @(x)[((exp(x(1))+exp(-x(1)))^2 - (exp(x(1))+exp(-x(1)))* ...
    (exp(x(1))-exp(-x(1))))/((exp(x(1))+exp(-x(1)))^2), 0; 0, 2];
Mkf = { hessf1, @(x) diag(diag(hessf1(x))), @(x) hessf1(x)};

%
Xka_hess = globalnewtonlike(f1,gradf1,Mkf{1},x_a,epsilon,maxiter,beta,gamma,a1,a2,p);
Xka_dia = globalnewtonlike(f1,gradf1,Mkf{2},x_a,epsilon,maxiter,beta,gamma,a1,a2,p);
Xka_approx = globalnewtonlike(f1,gradf1,Mkf{3}(x_a),x_a,epsilon,maxiter,beta,gamma,a1,a2,p);


%
Xkb_hess = globalnewtonlike(f1,gradf1,Mkf{1},x_b,epsilon,maxiter,beta,gamma,a1,a2,p);
Xkb_dia = globalnewtonlike(f1,gradf1,Mkf{2},x_b,epsilon,maxiter,beta,gamma,a1,a2,p);
Xkb_approx = globalnewtonlike(f1,gradf1,Mkf{3}(x_b),x_b,epsilon,maxiter,beta,gamma,a1,a2,p);

%
Xkc_hess = globalnewtonlike(f1,gradf1,Mkf{1},x_c,epsilon,maxiter,beta,gamma,a1,a2,p);
Xkc_dia = globalnewtonlike(f1,gradf1,Mkf{2},x_c,epsilon,maxiter,beta,gamma,a1,a2,p);
Xkc_approx = globalnewtonlike(f1,gradf1,Mkf{3}(x_c),x_c,epsilon,maxiter,beta,gamma,a1,a2,p);

figure(1);
clf; hold on
[x1, x2]  = meshgrid(-6:0.2:6, -6:0.2:6); 
fvals     = fv1(x1, x2);                      
hsurf     = surf(x1, x2, fvals);            
contours  = [0.001 0.05 0.1 0.5 1 3 5 10 20 30 100 1000];
[~,hcont] = contour(x1, x2, fvals, contours);
hsurf.FaceAlpha = 0.5;                       
hsurf.EdgeAlpha = 0.5;                        
hcont.LineWidth = 1.5;                       
xlabel('x1');                               
ylabel('x2'); 
title('GlobalNewtonLike for f1 and x_a')
plot3(Xka_hess(1,:),Xka_hess(2,:),0*Xka_hess(1,:),'k.-', 'LineWidth', 1.5, 'MarkerSize', 10.0)
plot3(Xka_dia(1,:),Xka_dia(2,:),0*Xka_dia(1,:),'r.-', 'LineWidth', 1.5, 'MarkerSize', 10.0)
plot3(Xka_approx(1,:),Xka_approx(2,:),0*Xka_approx(1,:),'g.-', 'LineWidth', 1.5, 'MarkerSize', 10.0)
saveas(gcf,'bild1.png')
hold on

figure(2)
clf; hold on
[x1, x2]  = meshgrid(-6:0.2:6, -6:0.2:6); 
fvals     = fv1(x1, x2);                      
hsurf     = surf(x1, x2, fvals);            
contours  = [0.001 0.05 0.1 0.5 1 3 5 10 20 30 100 1000];
[~,hcont] = contour(x1, x2, fvals, contours);
hsurf.FaceAlpha = 0.5;                       
hsurf.EdgeAlpha = 0.5;                        
hcont.LineWidth = 1.5;                       
xlabel('x1');                               
ylabel('x2'); 
title('GlobalNewtonLike for f1 and x_b')
plot3(Xkb_hess(1,:),Xkb_hess(2,:),0*Xkb_hess(1,:),'k.-', 'LineWidth', 1.5, 'MarkerSize', 10.0)
plot3(Xkb_dia(1,:),Xkb_dia(2,:),0*Xkb_dia(1,:),'r.-', 'LineWidth', 1.5, 'MarkerSize', 10.0)
plot3(Xkb_approx(1,:),Xkb_approx(2,:),0*Xkb_approx(1,:),'g.-', 'LineWidth', 1.5, 'MarkerSize', 10.0)
saveas(gcf,'bild2.png')
hold off

figure(3)
clf; hold on
[x1, x2]  = meshgrid(-6:0.2:6, -6:0.2:6); 
fvals     = fv1(x1, x2);                      
hsurf     = surf(x1, x2, fvals);            
contours  = [0.001 0.05 0.1 0.5 1 3 5 10 20 30 100 1000];
[~,hcont] = contour(x1, x2, fvals, contours);
hsurf.FaceAlpha = 0.5;                       
hsurf.EdgeAlpha = 0.5;                        
hcont.LineWidth = 1.5;                       
xlabel('x1');                               
ylabel('x2'); 
title('GlobalNewtonLike for f1 and x_c')
plot3(Xkc_hess(1,:),Xkc_hess(2,:),0*Xkc_hess(1,:),'k.-', 'LineWidth', 1.5, 'MarkerSize', 10.0)
plot3(Xkc_dia(1,:),Xkc_dia(2,:),0*Xkc_dia(1,:),'r.-', 'LineWidth', 1.5, 'MarkerSize', 10.0)
plot3(Xkc_approx(1,:),Xkc_approx(2,:),0*Xkc_approx(1,:),'g.-', 'LineWidth', 1.5, 'MarkerSize', 10.0)
saveas(gcf,'bild3.png')
hold off


%% Function 2
f2 = @(x) (x(1)^2 + x(2) - 11)^2 + (x(1) + x(2)^2 -7)^2;
fv2 = @(x1,x2) (x1.^2 + x2 -11).^2 + (x1 + x2.^2 - 7).^2;
gradf2 = @(x) [2*(x(1)^2 + x(2) - 11)*2*x(1) + 2*(x(1) + x(2)^2 -7); ...
    2*(x(1)^2 + x(2) - 11) + 2*(x(1) + x(2)^2 -7)*2*x(2)];
hessf2 = @(x) [12*x(1)^2 + 4*x(2)-42, 8*x(2); 4*x(1) + 4*x(2), 4*x(1) + 14*x(2)^2 - 26];

Mkf = { hessf1, @(x) diag(diag(hessf1(x))), @(x) hessf1(x)};
%
Xka_hess = globalnewtonlike(f2,gradf2,Mkf{1},x_a,epsilon,maxiter,beta,gamma,a1,a2,p);
Xka_dia = globalnewtonlike(f2,gradf2,Mkf{2},x_a,epsilon,maxiter,beta,gamma,a1,a2,p);
Xka_approx = globalnewtonlike(f2,gradf2,Mkf{3}(x_a),x_a,epsilon,maxiter,beta,gamma,a1,a2,p);


%
Xkb_hess = globalnewtonlike(f2,gradf2,Mkf{1},x_b,epsilon,maxiter,beta,gamma,a1,a2,p);
Xkb_dia = globalnewtonlike(f2,gradf2,Mkf{2},x_b,epsilon,maxiter,beta,gamma,a1,a2,p);
Xkb_approx = globalnewtonlike(f2,gradf2,Mkf{3}(x_b),x_b,epsilon,maxiter,beta,gamma,a1,a2,p);

%
Xkc_hess = globalnewtonlike(f2,gradf2,Mkf{1},x_c,epsilon,maxiter,beta,gamma,a1,a2,p);
Xkc_dia = globalnewtonlike(f2,gradf2,Mkf{2},x_c,epsilon,maxiter,beta,gamma,a1,a2,p);
Xkc_approx = globalnewtonlike(f2,gradf2,Mkf{3}(x_c),x_c,epsilon,maxiter,beta,gamma,a1,a2,p);

figure(4);
clf; hold on
[x1, x2]  = meshgrid(-6:0.2:6, -6:0.2:6); 
fvals     = fv2(x1, x2);                      
hsurf     = surf(x1, x2, fvals);            
contours  = [0.001 0.05 0.1 0.5 1 3 5 10 20 30 100 1000];
[~,hcont] = contour(x1, x2, fvals, contours);
hsurf.FaceAlpha = 0.5;                       
hsurf.EdgeAlpha = 0.5;                        
hcont.LineWidth = 1.5;                       
xlabel('x1');                         
ylabel('x2');
title('GlobalNewtonLike for f2 and x_a')
plot3(Xka_hess(1,:),Xka_hess(2,:),0*Xka_hess(1,:),'k.-', 'LineWidth', 1.5, 'MarkerSize', 10.0)
plot3(Xka_dia(1,:),Xka_dia(2,:),0*Xka_dia(1,:),'r.-', 'LineWidth', 1.5, 'MarkerSize', 10.0)
plot3(Xka_approx(1,:),Xka_approx(2,:),0*Xka_approx(1,:),'g.-', 'LineWidth', 1.5, 'MarkerSize', 10.0)
saveas(gcf,'bild4.png')
hold off

figure(5)
clf; hold on
[x1, x2]  = meshgrid(-6:0.2:6, -6:0.2:6); 
fvals     = fv2(x1, x2);                      
hsurf     = surf(x1, x2, fvals);            
contours  = [0.001 0.05 0.1 0.5 1 3 5 10 20 30 100 1000];
[~,hcont] = contour(x1, x2, fvals, contours);
hsurf.FaceAlpha = 0.5;                       
hsurf.EdgeAlpha = 0.5;                        
hcont.LineWidth = 1.5;                       
xlabel('x1');                               
ylabel('x2'); 
title('GlobalNewtonLike for f2 and x_b')
plot3(Xkb_hess(1,:),Xkb_hess(2,:),0*Xkb_hess(1,:),'k.-', 'LineWidth', 1.5, 'MarkerSize', 10.0)
plot3(Xkb_dia(1,:),Xkb_dia(2,:),0*Xkb_dia(1,:),'r.-', 'LineWidth', 1.5, 'MarkerSize', 10.0)
plot3(Xkb_approx(1,:),Xkb_approx(2,:),0*Xkb_approx(1,:),'g.-', 'LineWidth', 1.5, 'MarkerSize', 10.0)
saveas(gcf,'bild5.png')
hold off

figure(6)
clf; hold on
[x1, x2]  = meshgrid(-6:0.2:6, -6:0.2:6); 
fvals     = fv2(x1, x2);                      
hsurf     = surf(x1, x2, fvals);            
contours  = [0.001 0.05 0.1 0.5 1 3 5 10 20 30 100 1000];
[~,hcont] = contour(x1, x2, fvals, contours);
hsurf.FaceAlpha = 0.5;                       
hsurf.EdgeAlpha = 0.5;                        
hcont.LineWidth = 1.5;                       
xlabel('x1');                               
ylabel('x2'); 
title('GlobalNewtonLike for f2 and x_c')
plot3(Xkc_hess(1,:),Xkc_hess(2,:),0*Xkc_hess(1,:),'k.-', 'LineWidth', 1.5, 'MarkerSize', 10.0)
plot3(Xkc_dia(1,:),Xkc_dia(2,:),0*Xkc_dia(1,:),'r.-', 'LineWidth', 1.5, 'MarkerSize', 10.0)
plot3(Xkc_approx(1,:),Xkc_approx(2,:),0*Xkc_approx(1,:),'g.-', 'LineWidth', 1.5, 'MarkerSize', 10.0)
saveas(gcf,'bild6.png')
hold off


diary("off")

%% Results:
% One can observe that Hessian and the diagonal approximation have nearly
% the same results. For f1 they are identical because the hessian itself
% has zeros at its non-diagonal.
% (a) and (b) converge to their local minuma. As for f2 there are multiple
% local minima (picture). Method (c) on the other hand does not converge on
% muliple occasions, compare x_c for f2.