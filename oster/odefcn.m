function dydt = odefcn(y,p)
  dydt = zeros(2,1);
  dydt(1) = p(1)*y(1) - p(3)*y(1)*y(2);
  dydt(2) = -p(2)*y(2) + p(4)*y(1)*y(2);
end