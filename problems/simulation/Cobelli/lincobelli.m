function  [dx]=lincobelli(t,x,ut,u,par)

%u=zeros(size(0:0.01:1.27));
u=interp1(ut,u,t);

dx = zeros(2,1);    % a column vector
% dx(1) = -1/par(8)*(par(1))*x(1) - par(2)*x(2) + par(7)*u;
% dx(2) = 1/par(9)*(par(1)*x(1) - (par(5) + par(2))*x(2));
dx(1) = -(par(1))*x(1) + par(2)*x(2) + par(7)*u;
dx(2) = par(1)*x(1) - (par(5) + par(2))*x(2);

end