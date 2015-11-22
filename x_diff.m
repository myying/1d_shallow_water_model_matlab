function xo = x_diff( x, nx, order, type)
% type: 1=forward, -1=backward, 0=centered
% order: 2=2nd order; 4=4th order accurate
  if(order==2)
    j=2:nx-1;
    xo=x(j+1)-x(j-1);
  elseif(order==4)
    j=3:nx-2;
    xo=(4/3)*(x(j+1)-x(j-1))-(1/6)*(x(j+2)-x(j-2));
  end
end