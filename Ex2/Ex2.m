%Sheet 2 exercise 3


%% Part c
% Data
f = @(x) 100*(x(2) - x(1).^2).^2 + (1-x(1)).^2;
gradf = @(x,y) [200*(x(2) - x(1).^2).*(-2*x(1)) - 2*(1-x(1)); 200*(x(2) - x(1).^2)];
e = 10^-3;
%gamma = 10^-3;
%beta = 0.5;
maxit = 10000;
x0 = [ 1 -1.2 ; -0.5 1];
%x0  = [1; -0.5];
%x0 = [-1.2; 1];
%x0 = [-3; -3];
step = 0.2;
[Y1,Y2] = meshgrid(-3:step:3);
g = @(x,y) 100*(y - x.^2).^2 + (1-x).^2;
Z = g(Y1,Y2);

for j = 1:2
    figure(j)
    X = steepestdesc(f,gradf,x0(:,j),e,maxit);   
    surf(Y1,Y2,Z)
    hold on
    for i=1:length(X(1,:))
        plot3(X(1,i),X(2,i),f(X(:,i)),'X')
        hold on
    end
    hold off
    title(['x0: ', num2str(x0(1,j)),', ', num2str(x0(2,j))])
end
%% Part d
f = @(x) (2/3)*x.^3 + 0.5*x.^2;
gradf = @(x) 2*x.^2 + x;
x0 = [1 0.5 0.1];
%x0 = 1;
%x0 = 0.5;
%x0 = 0.1;
for j =1:3
    X = steepestdesc(f,gradf,x0(j),e,maxit);
    x = -5:step:5;
    figure(j+2)
    plot(x,f(x))
    hold on
    plot(X,f(X),'r*')
    hold off
    title(['x0: ', num2str(x0(j))])
end