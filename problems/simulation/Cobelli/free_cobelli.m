function [e,err]=free_cobelli(par)
 
%   loadfile = ['Results/fmindata.mat'];
%   load(loadfile)
 
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
par1=[k21,k12,Vm,Km,k02,c1,b1,V1,VE];
options = odeset('RelTol',1e-4,'AbsTol',[1e-4 1e-4]);

n=128;

%MODEL SET AS SYSTEM
    % [t,x0] = ode45(@cobelli1,[0:0.1:12.7],[x_i 0],options,par1);
     [t,x0] = ode45(@lincobelli,[0:0.1:12.7],[x_i 0],options,par1);

     k211 = par(1);
     k121 = par(2);
     k021 = par(3);
     
     par2=[k211,k121,Vm,Km,k021,c1,b1,V1,VE];
    % [t,x1] = ode45(@cobelli1,[0:0.1:12.7],[x_i 0],options,par2);
     [t,x1] = ode45(@lincobelli,[0:0.1:12.7],[x_i 0],options,par2);
    
    e=x0(1:n,1)-x1(1:n,1); % this is the simulation error

    err=sum(abs(e));  % Divide by number of points in error array

% e_save = [e_save e];
% err_save=[err_save err];
% par_save=[par_save par];

% fprintf(1,'Iteration %d: m=%.1f, c=%.3f d=%.1f --> Error = %5.3f\n',length(e_save),par(1),par(2),par(3),err);
% savefile = ['Results/fmindata.mat']; 

% savefile = ['Results/fmindata.mat'];
% save(savefile, 'e_save', 'err_save', 'par_save', '-mat') ;    