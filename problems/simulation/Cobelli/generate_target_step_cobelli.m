% generate cobelli step target

% define parameters
k21=3.033; 
k12=2.4;
Vm=378 ;
Km=2.88;
k02=1.5;
c1=1.0;
b1=1.0;
V1 = 79;
VE = 9;
par=[k21,k12,Vm,Km,k02,c1,b1,V1,VE];

% define parameters
params.simulation.ndata = 256;
params.simulation.time_span = .635*50;
params.simulation.samp_time=params.simulation.time_span/(params.simulation.ndata-1);
params.simulation.input_set=2;

pass=1;
ut = [0:params.simulation.samp_time:params.simulation.time_span];
t_span = [0 params.simulation.time_span];
simin=[ut', zeros(length(ut),1)];

%step function
simin(:,2) = [0; 100*ones(length(ut)-1,1)];
x_i=0;


options = simset('SrcWorkspace','current');

sim('lincobelli_mdl',t_span,options);
figure;    
plot(y_out);

dlmwrite('lincobelli_step_response.txt',y_out)
