function dxdt = lincob_str(t,x,ut,u,outstr,par)
%Bill LaCava 2013

u = interp1(ut,u,t); %interpolate data set (ut,u) at times t

outstr = regexprep(outstr,'u(1)','x(1)');
outstr = regexprep(outstr,'u(2)','x(2)');
outstr = regexprep(outstr,'u(3)','u');

dxdt = zeros(2,1);    % a column vector
eval(['dxdt(1) = ' outstr ';']);
dxdt(2) = (par(1)*x(1) - (par(5) + par(2))*x(2));

end