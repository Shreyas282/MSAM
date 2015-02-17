%% controller tuning

clear 

Sys = 'DC1';

for opt =0
for devcnt = 0
for rndcntr=1

% close all

Save=0;
% randcount=29;
% deviation=.75;
% p.wt_estimate=0;
    randcount = rndcntr;
    deviation = devcnt;
    p.wt_estimate=opt;
    continueontosecondphase=0;
    gamma= [1.0; 1.0];
    nom_gamma = gamma;
    %% number of iterations
    its = 3; % number of iterations for each perturbation case
    %%
    its2 = 2; % number of iterations for selected perturbation set
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

%  original parameters
% p.cons = {'Tk', 'Ti','Td';
%         4.0588,9.25,2.3125};
% tuned parameters
% p.cons = {'Tk', 'Ti','Td';
%         3.5294,16.800,4.200};

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
    
elseif strcmp(Sys,'DC2')
    p.simulation.sim_model='Sys2_MSA';
    p.simulation.time_span=400;
    p.cons = {'Tk', 'Ti','Td';
         0.88,28.01,5.45};
     
elseif strcmp(Sys,'Sys3')
    p.simulation.sim_model='Sys3_MSA';
    p.simulation.time_span=1000;
    p.cons = {'Tk', 'Ti','Td';
         0.64,50.40,15.56};

elseif strcmp(Sys,'Sys4')
    p.simulation.sim_model='Sys4_MSA';
    p.simulation.time_span=200;
    p.cons = {'Tk', 'Ti','Td';
          2.5568    0.7933    1.7947};
end
    
for k=1:length(p.cons(1,:))
    eval(['syms ' p.cons{1,k}])
end

%target
% exp(1) = 1.75; r1 = exp(1);
% exp(2) = 3; r2 = exp(2);
% p.mod_adapt.exp = exp;
% target = {['1/m*(-c*u(1)^' num2str(exp(1)) '-k*u(2)^' num2str(exp(2)) '+u(3))']};
% star_mod.eqn_sym = -c*dx1^r1-k*x1^r2+u_in;

exp(1) = 2; r1 = exp(1);
exp(2) = 3; r2 = exp(2);
p.mod_adapt.exp = exp;
%  dx1/dt=-3x1x3 -2x2x3 -3x3x3
% 
% dx2/dt=-3x1x2 +x1x3 -3x2x3
% 
% dx3/dt=3x1x2 +3x1x3 -x2x3
% a=3;
% b=2;
% c=3;
% 
% d=3;
% e=1;
% f=3;
% 
% g=3;
% h=3;
% i=1;

% target = {'Tk*u(1) + Tk/Ti*u(2) + Tk*Td*u(3)';};
% star_mod.eqn_sym = 1/m*(-c*dx1*abs(dx1)^(exp(1)-1)-k*x1*abs(x1)^(exp(2)-1)+u_in);
% target = {['1/m*(-c*u(1)*(abs(u(1)))^' num2str(exp(1)-1) '-k*u(2)*abs(u(2))^' num2str(exp(2)-1) '+u(3))']};
% target = {['1/m*(-c*u(1)*abs(u(1))^' num2str(exp(1)-1) '-k*u(2)*abs(u(2))^' num2str(exp(2)-1) '+u(3))']};
% star_mod.eqn_sym = (-c*(dx1./abs(dx1))*(abs(dx1))^(exp(1))-k*(x1./abs(x1))*(abs(x1))^exp(2)+u_in);

% simulation options
% p.simulation.sim_model='Sys1_MSA';
p.simulation.ndata=128;
% p.simulation.time_span=100;
p.simulation.partitioning=0;
p.mod_adapt.output_dydtheta=0;
% p.simulation.IC=[5 0 1];
% p.simulation.ICs = [.7 0 0;
%                     0 0.5 0;
%                     0 0 .89;
%                     .7 0.5 0;
%                     0 0.5 .89;
%                     0.7 0 .89;
%                     1 1 1];
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

continuept2 = 0; %run more iterations on part 2
if ~continuept2
    
    % nominal model
    nom_mod.eqn_sym = Tk*P + Ti*I + Td*D;
%     nom_mod.eqn_sym=D^(2653526402685845139/2305843009213693952)*Td + I^(62698626398682113501/73786976294838206464)*Ti + D^(552543512929762569/2305843009213693952)*P*Tk;
 
    %nom_mod.eqn_sym = candidate_models(nomcount);
    nom_mod = getTerms(nom_mod,'mod',p);
    num_terms = length(regexp([nom_mod.terms(:).type],'int'));
    
    % nom_mod.eqn_sym = nom_mod.t1.val+nom_mod.t2.val+nom_mod.t3.val;
    nom_mod.eqn_sym = GetEqnSym(nom_mod);
    nom_mod.eqn_str = GetEqnStr_sym(nom_mod,p.allvars);
%     nom_mod.eqn_str='u(3)^(2653526402685845139/2305843009213693952)*Td + u(2)^(62698626398682113501/73786976294838206464)*Ti + u(3)^(552543512929762569/2305843009213693952)*u(1)*Tk';
    disp(['nominal model: ' char(nom_mod.eqn_sym)]);
    
    %% get symbolic parameters
    
    
    % for count2 = 1:length([nom_mod.terms])
    %
    % %      nom_mod.terms(count2).val = sym(subs(nom_mod.terms(count2).val,{p.cons{1,:}},cons));
    % try
    %      nom_mod.terms(count2).val = simple(sym(nom_mod.terms(count2).val));
    % catch
    %     keyboard
    % end
    %      try
    %
    %          constantval(count2) = coeffs(nom_mod.terms(count2).val);
    %          constants(count2) = eval(['sym(''a' num2str(count2) ''');']);
    %          p.cons{1,end+1} = char(constants(count2));
    %          p.cons{2,end} = double(constantval(count2));
    %
    %     nom_mod.terms(count2).val = subs(nom_mod.terms(count2).val,constantval,constants);
    %
    %     catch
    %                          keyboard
    %     end
    %
    % end
    % nom_mod.eqn_sym = GetEqnSym(nom_mod);
    % nom_mod.eqn_str = GetEqnStr_sym(nom_mod,p.allvars);
    
    
    gamma= ones(num_terms,1);
    nom_gamma = gamma;
    
    n = length(p.intvars)^num_terms; % number of tweaks to try (hill climbs)
%    n=1;
    %% model adaptation settings
    p.mod_adapt.beta_start = .1*ones(num_terms,1);
    init_beta = repmat(p.mod_adapt.beta_start,1,its);
    p.mod_adapt.maxbeta = 0.3;
    p.mod_adapt.minbeta = .01;
    p.mod_adapt.betastep = .01;
    % betahist = [beta,zeros(length(beta),its)];
    p.mod_adapt.minDMC = .1;
    p.mod_adapt.maxDMC = 1;
    p.mod_adapt.adapt_beta = 0;    
    %mu settings
    p.mod_adapt.adapt_mu = 0;
    p.mod_adapt.mustart = repmat(0.15,num_terms,1); %starting confidence
    init_mu = repmat(p.mod_adapt.mustart,1,its);
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
    
    
    %% plotting options
    movie=0;
    plotinloop=1;
    plotinloop2=1;
    colors = ['r','b','g','m','c'];
    
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
     try
        tmp = RunCtrlModel(nom_mod.eqn_str,p);
        
        yhat = tmp;
        
    catch %#ok<*CTCH>
        disp('WARNING: Nominal model did not simulate succesfully.');
        yhat=ones(128,1)*rand*10000;
    end
    
    
    %% initial error, correlation and string distance
    nom_error = sum(abs(p.Y-yhat(:,1)));
    tmp = corrcoef(p.Y,yhat(:,1));
    nom_corr = tmp(1,2);
    %nom_corr2(nomcount) = max(xcorr(p.Y,y_mod(:,nomcount),'coeff'));
%     tmpstr1 = regexprep(char(star_mod.eqn_sym),p.cons(1,:),'c');
%     tmpstr2 = regexprep(char(nom_mod.eqn_sym),p.cons(1,:),'c');
%     nom_strdist = strdist(char(star_mod.eqn_sym),char(nom_mod.eqn_sym));
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
    
    % max_corr(1:its) = 0;
    max_corr = nom_corr;
    best_corr=nom_corr;
    best_gamma = [0,0];
    best_error = nom_error;
    sum_abs_error=nom_error;
    best_mod = nom_mod; % best model from all runs
    y_best = yhat;
    pass(1:n) = 0;
    p.mod_adapt.bestpass =1;
    p.mod_adapt.valleycontrol = 1;
    win_criteria = .9995;
    bite=1; % whether to 'bite' on iteration attempts, given simulation
    leggy = {};
    leg={};
    
    phi = zeros(p.simulation.ndata,num_terms,its,n);
    error = zeros(p.simulation.ndata,its);
    gamma=zeros(num_terms,its,n); %reset can_mod = nom_mod
    rho = zeros(num_terms,its,n);
    beta=[];
    p_t = zeros(num_terms,its,n);
    p_w = zeros(num_terms,its,n);
    delM = zeros(num_terms,its,n);
    mu=[];
    chosen=struct();
    pertnum=[];
    p.mod_adapt.MC_NLS_plot =0;
    if p.mod_adapt.MC_NLS_plot
        p.true_rho=exp;
        % option to plot known MC* delta error estimate comparison
        p.mod_adapt.output_dydtheta=1;
        tmp = RunPertModels_sym(star_mod,[0 0],p);
        p.mod_adapt.MCstar = tmp.MC;
        p.mod_adapt.dydtheta_star = tmp.dydtheta;
        
    end
    %% run hill climb
       for a = 1:n
        %reset variables for each hill climb
        
        can_mod = nom_mod;
        beta(:,:,a)=init_beta;
        
        mu(:,:,a) = init_mu;
        %     nom_mod.eqn_sym = GetEqnSym(nom_mod);
        %     can_mod.eqn_str = GetEqnStr_sym(can_mod);
        
        % choose unique set of perturbation types
        while pass(a) ~= 1
            try
                [chosen,can_mod,pertnum(:,a)] = ChoosePertType_sym(can_mod,p);
                models(:,a) = [can_mod;chosen];
            catch
                keyboard
            end
            if a==1
                pass(a)=1;
            end
            for q=1:length(pertnum(1,:))-1
                if pertnum(:,a) == pertnum(:,q)
                    pass(a)=0;
                    break
                else
                    pass(a)=1;
                end
            end
        end
        
        if pass(a)
            t_opts=find(strcmpi({can_mod.terms(:).type},'int'));
            %choose starting value of rho
            for j=1:length(t_opts)
                % if the perturbation variable is already part of the term
                tmp=regexp([char(can_mod.terms(t_opts(j)).val)],char(p.intvars(pertnum(j,a))),'once');
                if sum(tmp)>0
                    % if the perturbation variable is already part of the
                    % perturbed term, start rho at the exponent of that term
                    rho(j,:,a)=nom_gamma(j);
                else
                    % otherwise we are perturbing with a new variable, so rho
                    % is the size of the perturbation
                    try
                        rho(j,:,a) = beta(j,:,a);
                    catch
                        keyboard
                    end
                    
                end
            end
            % run gradient descent for chosen perturbation set
            %         if a==1
            if plotinloop
                h = figure;
            end
            %         end
            for x = 2:its
                %% --- model adaptation routine
                %% ----- compute sensitivity Phi
                
                %       if Produce_phi==1
                try
                    
                    [p_t(:,x,a),phi(:,:,x,a),error(:,x,a),sum_abs_error(x,a),max_corr(x,a),outputs(x,a),p,pass2(x,a)] = Sensitivity_compute(models(:,a),beta(:,x-1,a),gamma(:,x-1,a),rho(:,x-1,a),mu(:,x-1,a),p_t(:,x-1,a),p);
                    if p.mod_adapt.MC_NLS_plot == 1
                        p.dydtheta=outputs(x,a).dydtheta;
                        p.MC = outputs(x,a).MC;
                    end
                    %            delM = p_t;
                catch
                    disp('SensitivityCompute() failed');
                    keyboard
                    pass2(x,a)=0;
                end
                
				 % PCA
                
                 testmat = phi(:,:,x,a)'*phi(:,:,x,a);
                %phi_index(x,a) = sqrt(min(eig(testmat))/max(eig(testmat)));
                
                %           if wt_estimate
                %
                %               WT = [1,2]; % if w=1 wavelet='DerGauss'; w=2 wavelet='Sombrero'; if w=3 wavelet='Gauss'; if w=4 wavelet='Morlet';
                %               output_focus = [1];
                %               thresh_range=[1.2, 1.4, 1.6, 1.8];
                %
                %               try
                %               [p_w(:,x,a),mu_adapt(:,x,a)] = PARSIM(error(:,x,a),phi(:,:,x,a),WT,thresh_range,mu(:,x-1,a),rho(:,x-1,a));
                %               delM = p_w;
                %               catch
                %                   %keyboard
                %                   pass2(x,a)=0;
                %               end
                %           end
                %         gradient.models{nomcount} = models;
                %         gradient.p_t{nomcount} = p_t;
                %         gradient.phi{nomcount} = phi;
                %         gradient.error{nomcount} = error;
                %         gradient.sum_error{nomcount} = sum_abs_error;
                %         gradient.corr{nomcount} = max_corr;
                %         gradient.outputs{nomcount} = outputs;
                %         gradient.p{nomcount} = p;
                %         gradient.pass2{nomcount} = pass2;
                
                %         save(savefile_phi, 'gradient', '-mat');
                
                
                %       else
                %           load(savefile_phi)
                %
                %           p_t(:,x,a) = gradient.p_t{nomcount}(:,x,a);
                %           phi(:,:,x,a) = gradient.phi{nomcount}(:,:,x,a);
                %           error(:,x,a) = gradient.error{nomcount}(:,x,a);
                %           sum_abs_error(x,a) = gradient.sum_error{nomcount}(x,a);
                %           max_corr(x,a) = gradient.corr{nomcount}(x,a);
                %           outputs(x,a) = gradient.outputs{nomcount}(x,a);
                %           p = gradient.p{nomcount};
                %           pass2(x,a) = gradient.pass2{nomcount}(x,a);
                %       end
                
                %% INSERT SHAPE ANALYSIS CODE HERE
                %keyboard
                
                if pass2(x,a)
                    % set outputs
                    y0(:,:,x) = outputs(x,a).y0;
                    y_all(:,:,:,x)=outputs(x,a).y_all;
                    y_ave(:,:,x)=outputs(x,a).y_ave;
                    dMC(:,x,a) = outputs(x,a).dMC;
                    
                    %% valley jumping control 
                    if p.mod_adapt.valleycontrol==1
                        mu(:,x,a) = AvoidValleyJump(x,gamma(:,:,a),mu(:,x,a));
                    end
                    
                    
                    %% ----- adapt mu
                    if p.mod_adapt.adapt_mu
                        mu(:,x,a) = GetMu(p.Y,y_ave(:,:,x),p,p_t(:,x,a),beta(:,x-1,a));
%                     else
%                         mu(:,x,a) = mu(:,x-1,a);
                    end
                    
                   
                   %% ----- calculate gamma
%                     gamma(:,x,a) = mu(:,x,a).*p_t(:,x,a).*rho(:,x-1,a);
                    gamma(:,x,a) = mu(:,x,a).*p_t(:,x,a);
                    rho(:,x,a) = rho(:,x-1,a)+gamma(:,x,a);
                    
                    %         ----- update model
                    models(:,a) = UpdateModels_sym(models(:,a),gamma(:,x,a),p);
                    
                    %% ----- update beta size
                    if p.mod_adapt.adapt_beta
                        for j=1:length(dMC(:,x,a))
                            if abs(dMC(j,x,a)) < p.mod_adapt.minDMC
                                beta(j,x,a) = beta(j,x-1,a)+p.mod_adapt.betastep;
                            elseif abs(dMC(j)) > p.mod_adapt.maxDMC
                                beta(j,x,a) = beta(j,x-1,a)-p.mod_adapt.betastep;
                                %         elseif beta(j,end) > p.mod_adapt.init_beta(j)
                                %             beta(j) = beta(j) - p.mod_adapt.betastep;
                            end
                            if beta(j,x,a) > p.mod_adapt.maxbeta
                                beta(j,x,a) = p.mod_adapt.maxbeta;
                            elseif beta(j,x,a) < p.mod_adapt.minbeta
                                beta(j) = p.mod_adapt.minbeta;
                            end
                        end
                    else
                        beta(:,x,a) = .1;
                    end
                    
                    %% --- Plot Routine
                    if plotinloop==1
                        phiplot = reshape(phi(:,:,x,a),[p.simulation.ndata,num_terms]);
                        %             rhoplot = reshape(rho(:,:,a),size(rho,1),size(rho,3));
                        try
                            perttitle=RunPlotRoutine(h,a,x,its,y_ave,p,colors,...
                                sum_abs_error(1:x,a),max_corr(1:x,a),...
                                pertnum(:,a),p_t(:,x,a),phiplot,rho(:,:,a),error(:,:,a));
                        catch
                            keyboard
                        end
                    end
                    %     end
                    %% --- Update best model
                    if p.mod_adapt.bestpass
                        %update best model
                        %     if sum_abs_error(a,end)<best_error
                        if max_corr(x,a)/sum_abs_error(x,a) > best_corr/best_error
                            
                            best_mod=models(1,a);
                            
                            best_error=sum_abs_error(x,a);
                            best_corr = max_corr(x,a);
                            y_best = y_ave(:,1,x);
                            best_rho = rho(:,x,a);
                            best_rho_hist = rho;
                            best_models = models(:,a);
                            best_pertnum=pertnum(:,a);
                            best_mu = mu(:,x,a);
                            best_beta = beta(:,x,a);
                            endits = x;
                            if plotinloop==1
                                best_perttitle = perttitle;
                            else
                                best_perttitle='';
                            end
                            continueontosecondphase=1;
                        end
                        
                    end
                end
            end
            clearvars y0 y_all y_ave
        end
        close 
    end
    if continueontosecondphase
    %% plot winner from first round
    if plotinloop
        figure;
        plot(p.Y,'b'); hold on;
        plot(y_best,'r');
        title(['Best from initial run' best_perttitle]);
        l=legend('$y^*$','$\hat{y}$');
        set(l,'interpreter','latex','Fontsize',12);
    end
    
    
    %% Gradient descent for best perturbation case, starting from best model
    % keyboard
    
    gamma2(1:num_terms,1:its2)=0; %reset can_mod = nom_mod
    pertnum2 = best_pertnum;
    % gamhist2=[gamma,zeros(length(gamma),its2)];
    bite=1;
    rho2 = best_rho;
    mu2 = [best_mu best_mu];
    mu2 = repmat(p.mod_adapt.mustart,1,its2);
    beta2 = zeros(length(best_beta),its2);
    beta2(:,1) = best_beta;
    sum_abs_error2 = best_error;
    max_corr2 = best_corr;
    can_mod = best_models(1);
    p_t2 = zeros(num_terms,its);
    % delM2=p_t2;
    pertlist='';
    itsstart=2;
    end

elseif continuept2
    fprintf(['Continuing iterations...\n']);
    itsstart=its2;
    its2=2*its2;
end
if continueontosecondphase
for count=1:length(best_pertnum)
    pertlist=strcat(pertlist,',',char(p.intvars(best_pertnum(count))));
end
if plotinloop2
    h=figure;
end
pause(1.0);
for x = itsstart:its2
    pass3=1;
    
   % p.mod_adapt.algebra = 0;
    %% --- model adaptation routine
    try
        
        [p_t2(:,x),phi2(:,:,x),error2(:,x),sum_abs_error2(x),max_corr2(x),outputs,p,pass3] = Sensitivity_compute(best_models,beta2(:,x-1),gamma2(:,x-1),rho2(:,x-1),mu2(:,x-1),p_t2(:,x-1),p);
        if p.mod_adapt.MC_NLS_plot == 1
            p.dydtheta=outputs(x).dydtheta;
            p.MC = outputs(x).MC;
        end
        
    catch
%         keyboard
        pass3=0;
        
    end
    % try
    % [best_models,beta2(:,x),gamma2(:,x),rho2(:,x),mu2(:,x),p_t2(:,x),phi2(x,:,:),error2(:,x),sum_abs_error2(x),...
    %     max_corr2(x),outputs,p,pass3] = ModelAdaptationRoutine(best_models,...
    %     beta2(:,x-1),gamma2(:,x-1),rho2(:,x-1),mu2(:,x-1),p_t2(:,x-1),p);
    % catch
    %     keyboardf
    % end
    if pass3
        % set outputs
        y02(:,:,x) = outputs.y0;
        y_all2(:,:,:,x)=outputs.y_all;
        y_ave2(:,:,x)=outputs.y_ave;
        dMC2(:,x) = outputs.dMC;
        
        % valley jumping control - only works with adaptable mu
        if p.mod_adapt.valleycontrol==1
            mu2(:,x) = AvoidValleyJump(x,gamma2,mu2(:,x-1));
        end
        
        
        %% ----- adapt mu
        if p.mod_adapt.adapt_mu
            mu2(:,x) = GetMu(p.Y,y_ave2(:,:,x),p,p_t2(:,x),beta2(:,x-1));
        else
            mu2(:,x) = mu2(:,x-1);
        end
        
        %% ----- calculate gamma
        gamma2(:,x) = mu2(:,x).*p_t2(:,x).*rho2(:,x-1);
        rho2(:,x) = rho2(:,x-1)+gamma2(:,x);
        
        %         ----- update model
        best_models = UpdateModels_sym(best_models,gamma2(:,x),p);
        
        %% ----- update beta size
        beta2(:,x) = beta2(:,x-1);
        if p.mod_adapt.adapt_beta
            for j=1:length(dMC2(:,x))
                if abs(dMC2(j,x)) < p.mod_adapt.minDMC
                    beta2(j,x) = beta2(j,x)+p.mod_adapt.betastep;
                elseif abs(dMC(j)) > p.mod_adapt.maxDMC
                    beta2(j,x) = beta2(j,x)-p.mod_adapt.betastep;
                    %         elseif beta(j,end) > p.mod_adapt.init_beta(j)
                    %             beta(j) = beta(j) - p.mod_adapt.betastep;
                end
                if beta2(j,x) > p.mod_adapt.maxbeta
                    beta2(j,x) = p.mod_adapt.maxbeta;
                elseif beta2(j,x) < p.mod_adapt.minbeta
                    beta2(j) = p.mod_adapt.minbeta;
                end
            end
        else
            beta2(:,x) = .1;
        end
        
        % if pass3
        %         y02(:,:,x) = outputs.y0;
        %         y_all2(:,:,:,x)=outputs.y_all;
        %         y_ave2(:,:,x)=outputs.y_ave;
        %% --- plot routine
        if plotinloop2
            phiplot = reshape(phi2(:,:,x),[p.simulation.ndata,num_terms]);
            perttitle=RunPlotRoutine(h,1,x,its2,y_ave2,p,colors,...
                sum_abs_error2(1:x)',max_corr2(1:x)',...
                pertnum2,p_t2(:,x),phiplot,rho2,error2);
        end
        
%% --- Update best model
            if max_corr2(x)/sum_abs_error2(x) > best_corr/best_error
        %       if max_corr(end) > best_corr
                best_mod = best_models(1);
        
                best_error=sum_abs_error2(x);
                best_corr = max_corr2(x);
                y_best = y_ave2(:,1,end);
                best_rho = rho2(:,x);
        
            end

			if max_corr2(x)>win_criteria
        %         best_gamma = gamhist2(:,x);
                its_FTW = x;
                break;
            end
    end
end

if Save==1
    
    if p.wt_estimate==0 
        savefile = ['Control_Data/Controller_' p.simulation.sim_model '_NLS_randcount' num2str(randcount) '_mu' ...
            num2str(p.mod_adapt.mustart(1)*100) '_dev' num2str(deviation*100) '_its' num2str(its) '-' num2str(its2)]
    else
        savefile = ['Control_Data/Controller_' p.simulation.sim_model '_PARSIM_randcount' num2str(randcount) '_mu' ...
            num2str(p.mod_adapt.mustart(1)*100) '_dev' num2str(deviation*100) '_its' num2str(its) '-' num2str(its2)]
    end    
    save(savefile); 
    best_mod.eqn_sym
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
end
% clearvars -except rndcntr devcnt opt
end
% keyboard
end
end
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

