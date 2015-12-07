%% controller tuning

clear 
clc


 %% number of iterations
        its = 10; % number of iterations for each perturbation case
        its2 = 5; % number of iterations for selected perturbation set


%% saving options
        p.save=0;
        

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
        syms P I
        p.intvars = [P,I]; % internal variables of the system
        p.absintvars = [abs(P),abs(I)];
        p.mod_adapt.useabs =1;
        p.extvars = []; % external input variables of the system
        p.allvars = [p.intvars p.extvars];
        num_vars = length(p.allvars);
        p.simulation.input_set = 0;
        p.simulation.IC=0;

%% target controller output


Sys = 'NL_valve_2';
  
if strcmp(Sys,'NL_valve_2')



        p.simulation.time_span=100;    
        p.simulation.ndata=101;
        n = p.simulation.ndata;
        p.simulation.samp_time=p.simulation.time_span/(n-1);

        Cm=1;
        p.simulation.Cm=Cm;
    
        time_span=p.simulation.time_span;
        samp_time=p.simulation.samp_time;
        p.nl_order=4;
        par=[0.1,0.1];
        sim('NL_valve');   
        p.Y = y_tmp;
        u_target = u;
        t_sim=t;
        p.simulation.sim_model='NL_valve_MSAM_2';
    
        par_1=[0.1,0.1];
        p.cons = {'Tk', 'Ti' ;
                  par_1(1),par_1(2)};

 end

% p.savepath = '/data'; 

for k=1:length(p.cons(1,:))
    eval(['syms ' p.cons{1,k}])
end

        exp(1) = 2; r1 = exp(1);
        exp(2) = 3; r2 = exp(2);
        p.mod_adapt.exp = exp;
        p.simulation.partitioning=0;
        p.mod_adapt.output_dydtheta=0;

        % get target output
        limit1=1;
        limit2=72;
        limit3=1;
        limit4=128;
        thresh=0.2;
        p.mod_adapt.algebra = 0;

        p.continuept2 = 0; %run more iterations on part 2

    
%% nominal model
    p.nom_mod.eqn_sym = Tk*P + Ti*I;
    p.nom_mod.eqn_form = p.nom_mod.eqn_sym;
    p.nom_mod = getTerms(p.nom_mod,'mod',p);
    p.num_terms = length(regexp([p.nom_mod.terms(:).type],'int'));
    p.nom_mod.eqn_sym = GetEqnSym(p.nom_mod);
    p.nom_mod.eqn_str = GetEqnStr_sym(p.nom_mod,p.allvars,p);
    disp(['nominal model: ' char(p.nom_mod.eqn_sym)]);

    %% model adaptation settings
    p.mod_adapt.beta_start = .1*ones(p.num_terms,1);
    p.init_beta = repmat(p.mod_adapt.beta_start,1,its);
    p.mod_adapt.maxbeta = 0.3;
    p.mod_adapt.minbeta = .01;
    p.mod_adapt.betastep = .1;
    % betahist = [beta,zeros(length(beta),its)];
    p.mod_adapt.minDMC = .1;
    p.mod_adapt.maxDMC = 1;
    p.mod_adapt.adapt_beta = 0;    
    %mu settings
    p.mod_adapt.adapt_mu = 0;
    p.mod_adapt.mustart = repmat(0.1,p.num_terms,1); %starting confidence
    p.init_mu = repmat(p.mod_adapt.mustart,1,its);
    p.mod_adapt.agg = [15, .05]; %mu change aggressiveness
    p.mod_adapt.mu_min = 0.01;
    %% model complexity setting
    p.mod_adapt.struct_err_opt = 'error';
    p.mod_adapt.corr_trigger = 1;
    p.mod_adapt.err_trigger = 5;
    p.mod_adapt.mc_opt='pertsize'; % model complexity measure type
    p.mod_adapt.mc_opt2=0;
    p.mod_adapt.param_pert = 0.01;
    p.mod_adapt.param_pert_sides = 1; %1 for one sided, 2 for two-sided
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
    p.plotinloop2=1;
    p.mod_adapt.MC_NLS_plot =0;

    warning('off','all');
    % warning;

%% nominal model output
    [p.yhat, u_effort_initial] = RunCtrlModel(p.nom_mod.eqn_str,p);

   
    
    %% Program MSAM

   
    parallel = 0;
    
for opt =0
    for devcnt = 0
        for rndcntr = 0
            p.randcount = rndcntr;
            deviation = devcnt;
            p.wt_estimate=opt;
    
            %% Run MSAM
            tic;
            if parallel==1 
                    if matlabpool('size') == 0% checking to see if my pool is already open
                         matlabpool local;
                    end
                        out = MSAMp(p,its,its2);
        %               disp(['Total CPU time: ' num2str(cputime-t)]);
                        matlabpool close;
            else
                    out = MSAM(p,its,its2);
            end
            toc;
        end
    end
end
%% simulate final model
[y_best, u_effort_final] = RunCtrlModel_sym(out.best_mod,p);
figure; 
subplot(121)
plot(p.Y,'b'); hold on;
plot(p.yhat,'g');
plot(y_best,'r'); 
title('system response');
legend('target','initial model','final model');
subplot(122)
plot(u_target); hold on;
plot(u_effort_initial,'g'); 
plot(u_effort_final,'r'); 
legend('target controller','initial controller','final controller');
title('control effort');
% if p.save==1
% 
% mkdir(p.savepath);
% savefile = [p.savepath  '/' char(Sys) '_its' num2str(its) '-' num2str(its2)];  
% save(savefile); 
% 
% end
