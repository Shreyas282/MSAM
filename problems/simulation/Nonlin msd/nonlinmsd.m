function dxdt = nonlinmsd(t,x,ut,outstr,p)
% Bill LaCava 2013
% ode solver for the nonlinear mass spring damper.

u = interp1(ut,u,t); %interpolate data set (ut,u) at times t
for count=1:length([p.cons(1,:)])
        eval([p.cons{1,count} ' = ' num2str(p.cons{2,count}) ';']);
end
dxdt(1) = x(2);
eval(['dxdt(2) = ' outstr ';']);

end