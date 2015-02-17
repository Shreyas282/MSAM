function [pass,y_out,e,sume] = cobelli_error_mdl(outstr1,par,params)
% Bill LaCava 2013
% run candidate equation and calculate error versus correct cobelli.
% uses simulink instead of ode45
% outstr: candidate model
% par: cobelli parameters
% params: all the info you need
pass=1;
ut = [0:params.simulation.samp_time:params.simulation.time_span];
t_span = [0 params.simulation.time_span];
simin=[ut', zeros(length(ut),1)];

if params.simulation.input_set == 0 %free response
    simin(:,2) = zeros(length(ut),1);
    x_i=0.2;
elseif params.simulation.input_set == 1 %impulse
    simin(:,2) = [0;0;0;0;100;100;zeros(length(ut)-6,1)];
    x_i=0;
elseif params.simulation.input_set == 2 %step function
    simin(:,2) = [0; 100*ones(length(ut)-1,1)];
    x_i=0;
elseif params.simulation.input_set == 3 %random input
    simin(:,2) = wgn(length(ut),1,50);
    x_i=0;
end

% [t,x] = ode45(@(t,x) lincobelli(t,x,ut,u,par),t_span,[x_i 0],par);
% Y = x(:,1);
% clear t x
Y = params.target;
% options = odeset('RelTol',1e-4,'AbsTol',[1e-4 1e-4]);
options = simset('SrcWorkspace','current');

try
    sim('lincobelli_mdl',t_span,options);
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
    sume=0;
    y_out=e;
%     keyboard
    
end

end