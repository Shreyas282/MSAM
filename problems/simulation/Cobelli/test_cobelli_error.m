% test cobelli 

p.simulation.samp_time=.025;
p.simulation.time_span=.635*50;

p.simulation.input_set=2;

p.targetfile = 'lincobelli_step.txt';
tmp = textscan(fopen(p.targetfile,'rt'),'%f',1280);
p.target = tmp{:};

outstr='-par(1)*x(1) + par(2)*x(2) + par(7)*u';

x_i=0.2;

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

[pass,y,e,sume]=cobelli_error(outstr,par,p);

