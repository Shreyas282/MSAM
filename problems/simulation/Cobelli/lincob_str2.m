function dxdt = lincob_str2(t,x,ut,u,outstr,par)
%Bill LaCava 2013

u = interp1(ut,u,t); %interpolate data set (ut,u) at times t

dxdt = zeros(2,1);    % a column vector

dxdt(1) = par(7)*u-par(1)*x(1)+par(2)*x(2)
% dxdt(2) = (par(1)*x(1) - (par(5) + par(2))*x(2));
eval(['dxdt(2) = ' outstr ';']);

end