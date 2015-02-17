function [pass,y_out,e,sume] = cobelli_error(outstr,par,params)
% Bill LaCava 2013
% run candidate equation and calculate error versus correct cobelli.
% outstr: candidate model
% par: cobelli parameters
% params: all the info you need
warning('off','MATLAB:ode45:IntegrationTolNotMet')
pass=1;
ut = [0:params.simulation.samp_time:params.simulation.time_span];
t_span = ut;

if params.simulation.input_set == 0 %free response
    u = zeros(length(ut),1);
    x_i=0.2;
elseif params.simulation.input_set == 1 %impulse
    u = [0;0;0;0;100;100;zeros(length(ut)-6,1)];
    x_i=0;
elseif params.simulation.input_set == 2 %step function
    u = [0; 100*ones(length(ut)-1,1)];
    x_i=0;
elseif params.simulation.input_set == 3 %random input
    u = wgn(length(ut),1,50);
    x_i=0;
end
outstr = regexprep(outstr,'u\(1\)','x(1)');
outstr = regexprep(outstr,'u\(2\)','x(2)');
outstr = regexprep(outstr,'u\(3\)','u');

% [t,x] = ode45(@(t,x) lincobelli(t,x,ut,u,par),t_span,[x_i 0],par);
% Y = x(:,1);
% clear t x
Y = params.target;
% options = odeset('RelTol',1e-4,'AbsTol',[1e-4 1e-4]);

try
    [t,x] = ode45(@(t,x) lincob_str(t,x,ut,u,outstr,par),t_span,[x_i 0]);
    y_out=x(:,1);
catch
    pass=0;
    e=zeros(length(Y),1);
    sume=0;
    y_out=e;
%     keyboard
end
if pass==1 && length(y_out)==length(Y)
    try
    e = Y-y_out;
    sume=sum(abs(e));
    catch
%         keyboard
    end
else
    pass=0;
    e=zeros(length(Y),1);
    sume=e;
    y_out=e;
%     keyboard
    
end

end