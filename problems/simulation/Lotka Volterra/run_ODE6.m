% run lotka voltera
clear
x_in = [0:0.1:2];

time_span = 7;
samp_time=.025;
%figure;
options = simset('SrcWorkspace','current');
% for i=1:length(xs)
%     for j=1:length(ys)
%         x_i=xs(i);
%         y_i=ys(j);
%         try
%outstr='c1*u(1)-c2*u(1)*u(2)-c3*u(1)^2';
x=[];
ddy_all =[];
dy_all=[];
y_all=[];
x_all=[];
t_all=0;
for i=1:length(x_in)
    x = x_in(i);
   
    sim('ODE6.mdl');
    ddy_all=[ddy_all;ddy];
    dy_all=[dy_all;dy];
    y_all=[y_all;y];
    x_all=[x_all;repmat(x,length(ddy),1)];
    t_tmp = t_tmp+t_all(end);
    t_all = [t_all;t_tmp];
end
ode6 = [ddy_all,dy_all,y_all,x_all];

tmp = csaps(t_all(2:end),ddy_all);
sim_dy_all = trapz(ddy_all);
sim_y_all = trapz(y_all);

