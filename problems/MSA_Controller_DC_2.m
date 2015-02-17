%% controller tuning

clear 



 %% number of iterations
its = 3; % number of iterations for each perturbation case
its2 = 2; % number of iterations for selected perturbation set


%% saving options
p.save=0;
p.savepath = 'Control_Data/Controller_';

    continueontosecondphase=0;
    gamma= [1.0; 1.0];
    nom_gamma = gamma;
   
% Produce_eqns=1;
% savefile_models = 'Results\harm_eqns.mat';
% Produce_outs=1;
% savefile_outs = 'Results\harm_outputs.mat';
% Produce_phi=1;
% savefile_phi = 'Results\harm_phi_record.mat';
%% Model Setup
syms P I D 
p.intvars = [P, I D]; % internal variables of the system
% p.absintvars = [abs(x1),abs(x2),abs(x3)];
p.absintvars = [abs(P),abs(I),abs(D)];
p.mod_adapt.useabs = 0;
p.extvars = []; % external input variables of the system
p.allvars = [p.intvars p.extvars];
num_vars = length(p.allvars);
p.simulation.input_set = 0;
p.simulation.IC=0;

%% target controller output
Sys = 'DC1';
if strcmp(Sys,'DC1')
    
    
    p.simulation.time_span=15;
    p.simulation.R_filter = 10000; % 10K ohm
    p.simulation.C_filter = 30*10^-6; %30 microFarad
    p.simulation.num1 = [1, 0.71];
    p.simulation.den1= [4.4 3.7 1];
    p.simulation.ndata=128;
    n = p.simulation.ndata;
    p.simulation.sim_model='DC1_MSA';
    

    Tset = 2.3;
    %Zeta = 1.2;
    Zeta=1.2;
    %Mp = exp((-Zeta*pi)/(sqrt(1-Zeta^2)));
    Wn = 4/(Tset*Zeta);
    Targ_Spd = 2500; %Target motor speed in RPM
    Cm = 641.51; %Motor Speed conversion in RPM per volt
    Vref=Targ_Spd/Cm;
    p.simulation.Vref = Vref;
    sys = tf(Wn^2,[1 (2*Zeta*Wn) Wn^2]);
    t_s=p.simulation.time_span/(p.simulation.ndata-1);
    T = 0:t_s:p.simulation.time_span;

    p.Y = step(Vref*sys,T);
    %p.cons = {'Tk', 'Ti','Td';
    %    2.5496    0.7946    1.7958};
    p.cons = {'Tk', 'Ti','Td';
        5.5933    1.1417    4.5146};
    
% elseif strcmp(Sys,'DC2')
%     p.simulation.sim_model='Sys2_MSA';
%     p.simulation.time_span=400;
%     p.cons = {'Tk', 'Ti','Td';
%          0.88,28.01,5.45};
%      
% elseif strcmp(Sys,'Sys3')
%     p.simulation.sim_model='Sys3_MSA';
%     p.simulation.time_span=1000;
%     p.cons = {'Tk', 'Ti','Td';
%          0.64,50.40,15.56};
% 
% elseif strcmp(Sys,'Sys4')
%     p.simulation.sim_model='Sys4_MSA';
%     p.simulation.time_span=200;
%     p.cons = {'Tk', 'Ti','Td';
%           2.5568    0.7933    1.7947};
end
    
for k=1:length(p.cons(1,:))
    eval(['syms ' p.cons{1,k}])
end

exp(1) = 2; r1 = exp(1);
exp(2) = 3; r2 = exp(2);
p.mod_adapt.exp = exp;

p.simulation.ndata=128;
p.simulation.partitioning=0;
p.mod_adapt.output_dydtheta=0;

% get target output
limit1=1;
limit2=72;
limit3=1;
limit4=128;
thresh=0.2;

p.simulation.ndata=128;
n = p.simulation.ndata;
% p.simulation.time_span=4;
p.simulation.samp_time=p.simulation.time_span/(n-1);

% t=0:1:n-1;
% %load y_d1
% [num1,den1]=pade(10,20);
% sys1=tf(num1,den1);
% %sys2=tf([0.25],[1,1.0,0.25]);
% sys2=tf([1],[1,1]);
% sys3=series(sys1,sys2);
% %sys3=sys2;
% 
% p.Y = step(sys3,t);
% p = RunCtrlModel(target,p);

% samp_time=time_span/(n-1);
p.mod_adapt.algebra = 0;

p.continuept2 = 0; %run more iterations on part 2

    
%% nominal model
p.nom_mod.eqn_sym = Tk*P + Ti*I + Td*D;
p.nom_mod = getTerms(p.nom_mod,'mod',p);
p.num_terms = length(regexp([p.nom_mod.terms(:).type],'int'));
p.nom_mod.eqn_sym = GetEqnSym(p.nom_mod);
p.nom_mod.eqn_str = GetEqnStr_sym(p.nom_mod,p.allvars);
disp(['nominal model: ' char(p.nom_mod.eqn_sym)]);
%% model adaptation settings
p.mod_adapt.beta_start = .1*ones(p.num_terms,1);
p.init_beta = repmat(p.mod_adapt.beta_start,1,its);
p.mod_adapt.maxbeta = 0.3;
p.mod_adapt.minbeta = .01;
p.mod_adapt.betastep = .01;
% betahist = [beta,zeros(length(beta),its)];
p.mod_adapt.minDMC = .1;
p.mod_adapt.maxDMC = 1;
p.mod_adapt.adapt_beta = 0;    
%mu settings
p.mod_adapt.adapt_mu = 0;
p.mod_adapt.mustart = repmat(0.15,p.num_terms,1); %starting confidence
p.init_mu = repmat(p.mod_adapt.mustart,1,its);
p.mod_adapt.agg = [15, .05]; %mu change aggressiveness
p.mod_adapt.mu_min = 0.01;
%% model complexity setting
p.mod_adapt.struct_err_opt = 'error';
p.mod_adapt.corr_trigger = 1;
p.mod_adapt.err_trigger = 5;
p.mod_adapt.mc_opt='pardif'; % model complexity measure type
p.mod_adapt.mc_opt2=0;
% options for p.mc_opt:
% 'pardif' - dM = norm of difference in parameter sensitivity
% 'pertsize' - dM = beta (perturbation size)
% 'var' - dM = difference in variance given range of parameters
% 'varcorr' - dM = difference in correlation of variance given range of
%                   parameters
if strcmpi(p.mod_adapt.mc_opt,'pertsize')
    p.mod_adapt.adapt_beta=0;
end
p.mod_adapt.bestpass =1;
p.mod_adapt.valleycontrol = 1;
%% plotting options
p.plotinloop=1;
p.plotinloop2=0;
p.mod_adapt.MC_NLS_plot =0;
    
warning('off','all');
warning;
%% parameter error    

%     load('random3.mat')
%     %randcount = 3;
%     poodleydoo = 1;
%     for erw=((p.simulation.targetstate-1)*3+1):((p.simulation.targetstate-1)*3+1+2)
%         p.cons{2,erw} = p.cons{2,erw}*(1+deviation*(random(poodleydoo,randcount)));
%         poodleydoo=poodleydoo+1;
%     end
% %     p.cons = {'m', 'c', 'k';
% %         375*(1+deviation*(random(1,randcount))), 9800*(1+deviation*(random(2,randcount))), 130000*(1+deviation*(random(3,randcount)))};
%
%      try
%% nominal model output
p.yhat = RunCtrlModel(p.nom_mod.eqn_str,p);

%     catch %#ok<*CTCH>
%         disp('WARNING: Nominal model did not simulate succesfully.');
%         p.yhat=ones(128,1)*rand*10000;
%     end
    
    
    %% INITIAL MODEL PLOT
%     jbm = figure;
%     figure(jbm);
%     subplot(121)
%     plot(p.Y,'k'); hold on;
%     plot(yhat,'xk');
%     l=legend('$y^*(M^*,\Theta^*,t)$','Initial $\hat{y}(\hat{M},\hat{\Theta},t)$');
%     set(l,'interpreter','latex','fontsize',14);
%     xlabel('t','fontsize',14);
% % 
   

    
for opt =0
    for devcnt = 0
        for rndcntr = 0
            p.randcount = rndcntr;
            deviation = devcnt;
            p.wt_estimate=opt;
    
            %% Run MSAM
            out = MSAM(p,its,its2);
        end
    end
end
% end
%% Plot Best
% figure;
% plot(p.Y,'b'); hold on;
% plot(y_best(:,1,end),'r');
%set(gca,'fontsize',14);
% title('Best from gradient of winner from hill-climb');
% l=legend('$y^*$','$\hat{y}$');
% set(l,'interpreter','latex','Fontsize',12);
 %% FINAL MODEL PLOT
%  if rndcntr==1
%     figure(jbm);
%     subplot(122)
%     plot(p.Y,'k'); hold on;
%     plot(y_best,'xk');
%     l=legend('$y^*(M^*,\Theta^*,t)$','Final $\hat{y}(\hat{M},\Theta^*,t)$');
%     set(l,'interpreter','latex','fontsize',14);
%     xlabel('t','fontsize',14);
% %  end
% keyboard

% clearvars -except rndcntr devcnt opt

% keyboard


% keyboard
%end
%
% figure;
% subplot(121);
% surf(1:its2,1:p.simulation.ndata,phi2(:,:,1)');
% title('$\Phi _1$','Interpreter','Latex');
% xlabel('Iterations','Interpreter','Latex');
% ylabel('Time','Interpreter','Latex');
% subplot(122);
% surf(1:its2,1:p.simulation.ndata,phi2(:,:,2)');
% title('$\Phi _2$','Interpreter','Latex');
% xlabel('Iterations','Interpreter','Latex');
% ylabel('Time','Interpreter','Latex');
% end

