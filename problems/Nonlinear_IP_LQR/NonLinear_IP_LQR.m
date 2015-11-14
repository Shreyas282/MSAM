%% controller tuning

clear 
clc

 %% number of iterations
        its = 3; % number of iterations for each perturbation case
        its2 =25; % number of iterations for selected perturbation set
        p.gg =18;

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
     
        %syms P I D X
        syms x v theta w
        
        p.intvars = [x,v,theta,w]; % internal variables of the system
        p.absintvars = [abs(x),abs(v),abs(theta),abs(w)];
        p.mod_adapt.useabs = 0;
        p.extvars = []; % external input variables of the system
        p.allvars = [p.intvars p.extvars];
        num_vars = length(p.allvars);
        p.simulation.input_set = 0;
        p.simulation.IC=0;

%% target controller output


Sys = 'NonLinearPend_LQR';

    if strcmp(Sys,'NonLinearPend_LQR')
    
        p.simulation.M = 0.9;     % (kg) mass of cart
        p.simulation.m = 0.1;     % (kg) mass of pendulum
        p.simulation.b = 0;     % (N/m/s) friction of cart
        p.simulation.Lp = 0.47;   % (m) length of pendulum
        p.simulation.L = 0.235;   % (m) length of pendulum to center of gravity
        p.simulation.I = 0.0053;  % (kg-m^2) moment of inertia (pendulum)
        p.simulation.r = 0.023;   % (m) radius of pulley
        p.simulation.tau_m = 0.5; % (s) time constant of motor
        p.simulation.k_m = 17;    % (rad/s/V) gain of motor
        p.simulation.k_F = 9/pi;  % (rad/s/V) gain of feedback
        p.simulation.g = 9.81;    % (m/s^2) gravity
        p.simulation.Kc = 1;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Generating the Matrices and state space equation for pole      %
        % placement of IP                                                %  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         p.simulation.p = p.simulation.I*(p.simulation.M+p.simulation.m)+p.simulation.M*p.simulation.m*p.simulation.L^2; %denominator for the A and B matrices

            A = [0      1              0           0;
                0 -(p.simulation.I+p.simulation.m*p.simulation.L^2)*p.simulation.b/p.simulation.p  (p.simulation.m^2*p.simulation.g*p.simulation.L^2)/p.simulation.p   0;
                0      0              0           1;
                0 -(p.simulation.m*p.simulation.L*p.simulation.b)/p.simulation.p       p.simulation.m*p.simulation.g*p.simulation.L*(p.simulation.M+p.simulation.m)/p.simulation.p  0];
           
            B = [     0;
                      (p.simulation.I+p.simulation.m*p.simulation.L^2)/p.simulation.p;
                      0;
                      p.simulation.m*p.simulation.L/p.simulation.p
                ];
            
            C = [1 0 0 0;
                0 0 1 0];
            
            D1 = [0;0];

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        p.simulation.time_span=10;
        p.simulation.ndata=101;
        p.simulation.samp_time=p.simulation.time_span/(p.simulation.ndata-1);
       %p.simulation.sim_model='Pend_MSAM_PP'; %Pend_MSAM_PP';
        p.simulation.sim_model='Pend_MSAM_PP_noObserver';
 %% System TF for use in Runctrlmodel   
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         p.simulation.A=A;
         p.simulation.B=B;
         p.simulation.C=C;
         p.simulation.D=D1;
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %Generating target response using Pole placement of linear model and     %
 %an Observer                                                             % 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Pole Placement    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Pole placement gain
            p1 = -1;
            p2 = -2;
            p3 = -4.7379;
            p4 = -4.7378;
% Observer gain 
            l1=-1;
            l2=-1.1;
            l3=-1.2;
            l4=-1.3;
            
            K = place(p.simulation.A,p.simulation.B,[p1 p2 p3 p4]);
            L=place(p.simulation.A',p.simulation.C',[l1,l2,l3 l4])';
            
 %%%%%%%%%%%%%%%%%%%%%%
 %        LQR         %
 %%%%%%%%%%%%%%%%%%%%%%
%             Q = C'*C;
%             R = 1;
%             K = lqr(A,B,Q,R);
           
  
p.cons = {'Tk', 'Ti','Td','Te';K(1),K(2),K(3),K(4)};

    for k=1:length(p.cons(1,:))
        eval(['syms ' p.cons{1,k}])
    end

    %%%%%%%%%%%% Observer %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            p.simulation.A_hat=p.simulation.A;
            p.simulation.B_hat=[p.simulation.B,L];
            p.simulation.C_hat=[1,0,0,0;0,1,0,0;0,0,1,0;0,0,0,1];
            p.simulation.D_hat=[0 0 0;0 0 0;0 0 0;0 0 0];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
    
            p.nom_mod.eqn_sym = Tk*x+Ti*v+Td*theta+Te*w;
            p.nom_mod.eqn_form = p.nom_mod.eqn_sym;
            p.nom_mod = getTerms(p.nom_mod,'mod',p);
            p.num_terms = length(regexp([p.nom_mod.terms(:).type],'int'));
            p.nom_mod.eqn_sym = GetEqnSym(p.nom_mod);
            p.nom_mod.eqn_str = GetEqnStr_sym(p.nom_mod,p.allvars);
            disp(['Target model: ' char(p.nom_mod.eqn_sym)]);

            y_tmp = RunCtrlModel(p.nom_mod.eqn_str,p);
            p.Y=y_tmp(:,2); %%%%%%%%%%%%% this is the angle of the IP
            a=p.Y; % this is the target. 
            
           % save(['Y_' num2str(p.simulation.time_span) 'sec'],a);
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %                          Non Linear Model 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
  
       % parameters for non linear model 
        p.simulation.k1=p.simulation.M+p.simulation.m;
        p.simulation.k2=p.simulation.b; 
        p.simulation.k3=p.simulation.m*p.simulation.L;
        p.simulation.k4=(p.simulation.I+p.simulation.m*p.simulation.L^2); 
        p.simulation.k5=p.simulation.m*p.simulation.g*p.simulation.L;
        

        %Change model to simulate non linear IP
         p.simulation.sim_model='NonLinear_IPend_PP_Michigan_MSAM';

         y_tmp = RunCtrlModel(p.nom_mod.eqn_str,p); % this is taking too long to simulate
         p.Y2=y_tmp;
        
        
    
        plot(y_tmp,'b');hold on; plot(p.Y,'r*');
              
              
    end

% p.folderName = ['data/' char(Sys)];    

% p.savepath = ['data/' char(Sys)] ;
 p.savepath = ['data/'];
    
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
    p.nom_mod.eqn_sym = Tk*x+Ti*v+Td*theta+Te*w ;
    p.nom_mod.eqn_form = p.nom_mod.eqn_sym;
    p.nom_mod = getTerms(p.nom_mod,'mod',p);
    p.num_terms = length(regexp([p.nom_mod.terms(:).type],'int'));
    p.nom_mod.eqn_sym = GetEqnSym(p.nom_mod);
    p.nom_mod.eqn_str = GetEqnStr_sym(p.nom_mod,p.allvars);
    disp(['nominal model: ' char(p.nom_mod.eqn_sym)]);

    %% model adaptation settings
    p.mod_adapt.beta_start = .2*ones(p.num_terms,1);
    p.init_beta = repmat(p.mod_adapt.beta_start,1,its);
    p.mod_adapt.maxbeta = 0.3;
    p.mod_adapt.minbeta = .01;
    p.mod_adapt.betastep = .2;
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
    p.plotinloop2=1;
    p.mod_adapt.MC_NLS_plot =0;

    warning('off','all');
    % warning;

%% nominal model output
%     set_param('NonlinearPendulum','AlgebraicLoopSolver','LineSearch');
    p.yhat = RunCtrlModel(p.nom_mod.eqn_str,p);

%    p.Y = zeros(size(p.yhat));
    
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
                        out = MSAMp_partial(p,its,its2);
        %               disp(['Total CPU time: ' num2str(cputime-t)]);
                        matlabpool close;
            else
                    out = MSAM_partial(p,its,its2);
            end
            toc;
        end
    end
end

if p.save==1

% mkdir(p.savepath);
savefile = [p.savepath  '/' char(Sys) '_its' num2str(its) '-' num2str(its2)];  
save(savefile); 

end
