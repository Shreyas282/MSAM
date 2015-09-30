%% VIV force equation

clear 



 %% number of iterations
its = 3; % number of iterations for each perturbation case
its2 = 10; % number of iterations for selected perturbation set


%% saving options
p.save=1;


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
syms QDDOT QDOT Q
p.intvars = [QDDOT QDOT Q]; % internal variables of the system
% p.absintvars = [abs(x1),abs(x2),abs(x3)];
p.absintvars = [abs(QDDOT),abs(QDOT),abs(Q)];
p.mod_adapt.useabs = 1;
p.extvars = []; % external input variables of the system
p.allvars = [p.intvars p.extvars];
num_vars = length(p.allvars);

%% load target
viv = load('problems\VIV\VIV_data.mat');
for s=19
    p.savepath = ['D:\MSAM\VIV\VIV_' num2str(s)];
set = s;

tmp_cons = par_estimate_VIV(set);
% A_val = 12;
% eps_val = 0.0138;
p.cons = {'A', 'eps', 'St', 'D','U'; 
           12,  0.3, viv.St, viv.D, viv.U(set)};
% p.cons = {'B', 'C', 'E'; 
% %           -1.1584 , -2.7181, -217.4027};
%       2.2732    -0.4136    273.9277};
% make parameters symbolic
for k=1:length(p.cons(1,:))
    eval(['syms ' p.cons{1,k} ';']);
end
% p.Y = -viv.A/viv.D*viv.ddy_50s(set,:)';
% p.qddot = viv.ddforce_50s(set,:)';
% p.qdot = viv.dforce_50s(set,:)';
% p.q = viv.force_50s(set,:)';
train = 1:ceil(size(viv.ddy_50s,2)/2);

p.Y = -p.cons{1,1}/viv.D*reshape(viv.ddy_50s(set,train)',size(viv.ddy_50s(set,train)',1)*size(viv.ddy_50s(set,train)',2),1);
p.QDDOT = reshape(viv.ddforce_50s(set,train)',size(viv.ddforce_50s(set,train)',1)*size(viv.ddforce_50s(set,train)',2),1);
p.QDOT = reshape(viv.dforce_50s(set,train)',size(viv.dforce_50s(set,train)',1)*size(viv.dforce_50s(set,train)',2),1);
p.Q = reshape(viv.force_50s(set,train)',size(viv.force_50s(set,train)',1)*size(viv.force_50s(set,train)',2),1);

p.simulation.ndata = length(p.Q);

    
% p.mod_adapt.exp = exp;

% p.simulation.ndata=128;
p.simulation.partitioning=0;
p.mod_adapt.output_dydtheta=0;

% p.simulation.ndata=128;
% n = p.simulation.ndata;
% p.simulation.time_span=4;
% p.simulation.samp_time=p.simulation.time_span/(n-1);

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
p.mod_adapt.algebra = 1;

p.continuept2 = 0; %run more iterations on part 2

    
%% nominal model
p.nom_mod.eqn_sym = (QDDOT + eps*(2*pi*St*U/D)*(Q^2-1)*QDOT...
                    + (2*pi*St*U/D)^2*Q);
% p.nom_mod.eqn_sym = B*QDDOT + C*QDOT+E*Q;
p.nom_mod = getTerms(p.nom_mod,'mod',p);
p.num_terms = length(regexp([p.nom_mod.terms(:).type],'int'));
p.nom_mod.eqn_sym = GetEqnSym(p.nom_mod);
p.nom_mod.eqn_str = GetEqnStr_sym(p.nom_mod,p.allvars);
p.nom_mod.eqn_form = GetEqnForm(p.nom_mod);
disp(['nominal model: ' char(p.nom_mod.eqn_form)]);
%% model adaptation settings
p.mod_adapt.beta_start = .01*ones(p.num_terms,1);
p.init_beta = repmat(p.mod_adapt.beta_start,1,its);
p.mod_adapt.maxbeta = 0.1;
p.mod_adapt.minbeta = .001;
p.mod_adapt.betastep = .001;
% betahist = [beta,zeros(length(beta),its)];
p.mod_adapt.minDMC = .1;
p.mod_adapt.maxDMC = 1;
p.mod_adapt.adapt_beta = 0;    
%mu settings
p.mod_adapt.adapt_mu = 0;
p.mod_adapt.mustart = repmat(0.5,p.num_terms,1); %starting confidence
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
p.plotinloop=0;
p.plotinloop2=0;
p.mod_adapt.MC_NLS_plot =0;
    
warning('off','all');
% warning;
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
p.yhat = EvalModel(p.nom_mod.eqn_sym,p);

%     catch %#ok<*CTCH>
%         disp('WARNING: Nominal model did not simulate succesfully.');
%         p.yhat=ones(128,1)*rand*10000;
%     end
    
    
    %% INITIAL MODEL PLOT
%     jbm = figure;
%     figure(jbm);
%     subplot(121)
%     plot(p.Y,'k'); hold on;
%     plot(p.yhat,'xk');
%     l=legend('$y^*(M^*,\Theta^*,t)$','Initial $\hat{y}(\hat{M},\hat{\Theta},t)$');
%     set(l,'interpreter','latex','fontsize',14);
%     xlabel('t','fontsize',14);
% % 
   
parallel = 1;
    
for opt =0
    for devcnt = 0
        for rndcntr = 0
            p.randcount = rndcntr;
            deviation = devcnt;
            p.wt_estimate=opt;
    
            %% Run MSAM
            tic;
            if parallel 
                if matlabpool('size') == 0% checking to see if my pool is already open
                 matlabpool local;
                end
                out = MSAMp(p,its,its2);
%                 disp(['Total CPU time: ' num2str(cputime-t)]);
%                 matlabpool close;
            else
                out = MSAM(p,its,its2);
            end
            toc;
        end
    end
end
end
% end
%% tune parameters of final model
cons_tuned = par_estimate_VIV(19,{char(out.best_mod.eqn_sym)});
tmp = p.cons;
p.cons(1:2,1:2) = cons_tuned;
ybest_tuned = EvalModel(out.best_mod.eqn_sym,p);
p.cons = tmp;
%% Plot validation result
test = (train(end)+1):(size(viv.ddy_50s,2));
p.Y = -p.cons{1,1}/viv.D*reshape(viv.ddy_50s(set,test)',size(viv.ddy_50s(set,test)',1)*size(viv.ddy_50s(set,test)',2),1);
p.QDDOT = reshape(viv.ddforce_50s(set,test)',size(viv.ddforce_50s(set,test)',1)*size(viv.ddforce_50s(set,test)',2),1);
p.QDOT = reshape(viv.dforce_50s(set,test)',size(viv.dforce_50s(set,test)',1)*size(viv.dforce_50s(set,test)',2),1);
p.Q = reshape(viv.force_50s(set,test)',size(viv.force_50s(set,test)',1)*size(viv.force_50s(set,test)',2),1);

ybest_test = EvalModel(out.best_mod.eqn_sym,p);
ynom_test = EvalModel(out.nom_mod.eqn_sym,p);

% test parametrically tuned model on validation set
tmp = p.cons;
p.cons = cons_tuned;
ybest_tuned_test = EvalModel(out.best_mod.eqn_sym,p);
p.cons = tmp;

ybest_tuned_test = EvalModel(out.best_mod.eqn_sym,p);
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

