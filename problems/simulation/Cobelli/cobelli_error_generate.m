% This is the random error generator for the cobelli model
clear all

max_icount = 100;

e_save=[];
err_save=[];
par_save=[];

% savefile = ['Results/fmindata.mat'];
% save(savefile, 'e_save', 'err_save', 'par_save', '-mat')

load('cobelli_start.mat');
for icount = 1:max_icount
    icount
par = cobelli_start(:,icount);
[e(:,icount),err(icount)] = free_cobelli(par);
end

fname = ['Results/lin_lin_cobelli.mat'];
save(fname, 'e', 'err', '-mat');