function out = MSAMp(p,its,its2)
%Model Structure Adaptation Method
% William La Cava 2015
try
out=[];
disp('\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\');
disp(' Running MSAMp version 1.0.0-wgl');
disp('///////////////////////////////');

disp(['Number of round robin iterations: ' num2str(its)]);
disp(['Number of final adaptation iterations: ' num2str(its2)]);
digits(6);
if ~p.continuept2    
     colors = ['r','b','g','m','c'];
    %% initial error, correlation and string distance
    num_models = length(p.intvars)^p.num_terms; % number of models to try (hill climbs)
%     num_models=8;
%pert_index = zeros(length(p.intvars),p.num_terms);
    a = repmat({1:length(p.intvars)},1,p.num_terms);
    pert_index = allcomb(a{:});
    pert_index = shuffle(pert_index,1);
    nom_error = sum(abs(p.Y-p.yhat(:,1)));
    tmp = corrcoef(p.Y,p.yhat(:,1));
    nom_corr = tmp(1,2);
%     nom_corr = R2(p.Y,p.yhat(:,1));
    
    sum_abs_error=zeros(num_models,its); 
    max_corr=zeros(num_models,its); 
    best_corr_p  = nom_corr*ones(1,num_models);
    best_error_p = nom_error*ones(1,num_models);
    best_models(1:num_models) = p.nom_mod;
    best_mod = p.nom_mod; % best model from all runs
    best_rho = zeros(p.num_terms,1);
    best_rho_p = zeros(p.num_terms,num_models);
    best_mu_p = repmat(p.mod_adapt.mustart,1,num_models);
    best_beta_p = zeros(p.num_terms,num_models);
%     best_gamma = [0,0];
    y_best = p.yhat;
    y_best_p = repmat(p.yhat,1,num_models);
    pass(1:num_models) = 0;
    p.mod_adapt.bestpass =1;
    p.mod_adapt.valleycontrol = 1;
    win_criteria = .9995;
    best_models_p = repmat(p.nom_mod,p.num_terms+1,num_models);   
    nom_gamma = ones(p.num_terms,1);
    phi = zeros(num_models,p.simulation.ndata,p.num_terms,its);
   
    gamma=zeros(num_models,p.num_terms,its); %reset can_mod = nom_mod
    rho = zeros(num_models,p.num_terms,its);
    for i=1:num_models, beta(i,:,:)= p.init_beta; end
    p_t = zeros(num_models,p.num_terms,its);
    p_w = zeros(num_models,p.num_terms,its);
    delM = zeros(num_models,p.num_terms,its);
    mu=[];
%     best_perttitle_p(1:a) = '';
    chosen=struct();
    outputs= struct('y0',cell(num_models,its),'y_all',[],'y_ave',[],'dMC',[]);
   % pertnum=zeros(p.num_terms);
%    models = [];
    %% run hill climb
   
%    n=5;
% matlabpool;
warning('off','all');
parfor a = 1:num_models
    warning('off','all');

    %reset variables for each hill climb
    
    can_mod = p.nom_mod;
     
    %     nom_mod.eqn_sym = GetEqnSym(nom_mod);
    %     can_mod.eqn_str = GetEqnStr_sym(can_mod);

    % choose unique set of perturbation types
%         while pass(a) ~= 1
            [chosen,can_mod] = ChoosePertType(can_mod,p,pert_index(a,:));
            models = [can_mod;chosen];
    disp(['Running model adaptation ' num2str(a) ' of ' ...
           num2str(num_models) ': y = ' char(models(1).eqn_form)]);
%             if a==1
%                 pass(a)=1;
%             end
%             for q=1:length(pertnum(1,:))-1
%                 if pertnum(:,a) == pertnum(:,q)
%                     pass(a)=0;
%                     break
%                 else
%                     pass(a)=1;
%                 end
%             end
%         end

%         if pass(a)
        t_opts=find(strcmpi({can_mod.terms(:).type},'int'));
        %choose starting value of rho
        tmp_beta = beta(a,:,:);
        tmp_rho = rho(a,:,:);
        pert_index_a = pert_index(a,:);
        for j=1:length(t_opts)
            % if the perturbation variable is already part of the term
            tmp=regexp([char(can_mod.terms(t_opts(j)).val)],char(p.intvars(pert_index_a(j)')),'once');
            if sum(tmp)>0
                % if the perturbation variable is already part of the
                % perturbed term, start rho at the exponent of that term
                tmp_rho(j) = nom_gamma(j);
                %rho(a,j,:)=nom_gamma(j);
            else
                % otherwise we are perturbing with a new variable, so rho
                % is the size of the perturbation
                tmp_rho(:,j,:) = tmp_beta(:,j,:);
                
            end
        end
        rho(a,:,:) = tmp_rho;
        % run gradient descent for chosen perturbation set
        %         if a==1
        if p.plotinloop
            h = figure;
        end
        %         end
        p_t_a = p_t(a,:,:); 
        p_t_a = reshape(p_t_a,size(p_t_a,2),size(p_t_a,3));
        
        phi_a = phi(a,:,:,:);
        phi_a = reshape(phi_a,size(phi_a,2),size(phi_a,3),size(phi_a,4));
%         error_a = error(a,:,:);
        sum_abs_error_a = [];%sum_abs_error(a,:);
        max_corr_a = []; %max_corr(a,:);
       outputs_a = struct('y0',cell(1,its),'y_all',[],'y_ave',[],'dMC',[]);
%        outputs_a(1:its).y_all= zeros(p.simulation.ndata,1);
%        outputs_a(1:its).y_ave = zeros(p.simulation.ndata,1);
%        outputs_a(1:its).dMC = zeros(p.simulation.ndata,1);
       error = zeros(p.simulation.ndata,its);
       beta_a = p.init_beta;
       mu_a = p.init_mu;
%         beta_a = reshape(beta(a,:,:),size(beta,2),size(beta,3));
        gamma_a = gamma(a,:,:);
        gamma_a = reshape(gamma_a,size(gamma_a,2),size(gamma_a,3));
        rho_a = rho(a,:,:);
        rho_a = reshape(rho_a,size(rho_a,2),size(rho_a,3));
%         mu_a = reshape(mu(a,:,:),size(mu,2),size(mu,3));
%         best_corr_tmp = nom_corr;
%         best_erro_tmp = nom_error;
        dMC_a = zeros(p.num_terms,its);
        model_disp = 0;
        for x = 2:its
            
            %% --- model adaptation routine
            %% ----- compute sensitivity Phi

            %       if Produce_phi==1
%                 try
                
                [p_t_a(:,x),phi_a(:,:,x),error(:,x),sum_abs_error_a(x),max_corr_a(x),outputs_a(x),pass2] = Sensitivity_compute(models,beta_a(:,x-1),gamma_a(:,x-1),rho_a(:,x-1),p_t_a(:,x-1),p);
%                 if p.mod_adapt.MC_NLS_plot == 1
%                     p.dydtheta=outputs(a,x).dydtheta;
%                     p.MC = outputs(a,x).MC;
%                 end
                %            delM = p_t;
%                 catch
%                     disp('SensitivityCompute() failed');
%                     %keyboard
%                     pass2(x,a)=0;
%                 end

             % PCA

             testmat = phi_a(:,:,x)'*phi_a(:,:,x);
          

            if pass2
                % set outputs
               y0 = outputs_a(x).y0;
%                y_all=outputs_a(x).y_all;
               y_ave=outputs_a(x).y_ave;
                dMC_a(:,x) = outputs_a(x).dMC;
%% --- Update best model
                if p.mod_adapt.bestpass
                    %update best model
                    %     if sum_abs_error(a,end)<best_error
                    if max_corr_a(x)/sum_abs_error_a(x) > best_corr_p(a)/best_error_p(a)

%                         best_mod_p(a)=models(a,1);
%                         if ~model_disp
%                             disp(['best model updated: y = ' char(models(1).eqn_form)]); 
%                             model_disp=true; 
%                         end
                        best_error_p(a)=sum_abs_error_a(x);
                        best_corr_p(a) = max_corr_a(x);
                        y_best_p(:,a) = outputs_a(x).y0(:,1);
                        best_rho_p(:,a) = rho_a(:,x-1);
%                         best_rho_hist = rho_a;
                        best_models_p(:,a) = models;
                        best_pertnum_p(a,:) =pert_index(a,:);
                        best_mu_p(:,a) = mu_a(:,x-1);
                        best_beta_p(:,a) = beta_a(:,x);
%                         endits = x;
%                         if p.plotinloop==1
%                             best_perttitle_p(a) = perttitle;
%                         else
%                             best_perttitle_p(a)='';
%                         end
%                         continueontosecondphase(a)=1;
                    end

                end
                %% valley jumping control 
                if p.mod_adapt.valleycontrol==1
                    mu_a(:,x) = AvoidValleyJump(x,gamma_a(:,:),mu_a(:,x));
                end


                %% ----- adapt mu
                if p.mod_adapt.adapt_mu
                    mu_a(:,x) = GetMu(p.Y,y_ave,p,p_t_a(:,x),beta_a(:,x-1));
%                     else
%                         mu(:,x,a) = mu(:,x-1,a);
                end


               %% ----- calculate gamma
%                     gamma(:,x,a) = mu(:,x,a).*p_t(:,x,a).*rho(:,x-1,a);
                gamma_a(:,x) = mu_a(:,x).*p_t_a(:,x);
                rho_a(:,x) = rho_a(:,x-1)+gamma_a(:,x);

                %         ----- update model
                models = UpdateModels_sym(models,gamma_a(:,x),p);

                %% ----- update beta size
                if p.mod_adapt.adapt_beta
                    for j=1:length(dMC_a(:,x))
                        if abs(dMC_a(j,x)) < p.mod_adapt.minDMC
                            beta_a(j,x) = beta_a(j,x-1)+p.mod_adapt.betastep;
                        elseif abs(dMC_a(j,x)) > p.mod_adapt.maxDMC
                            beta_a(j,x) = beta_a(j,x-1)-p.mod_adapt.betastep;
                            %         elseif beta(j,end) > p.mod_adapt.p.init_beta(j)
                            %             beta(j) = beta(j) - p.mod_adapt.betastep;
                        end
                        if beta_a(j,x) > p.mod_adapt.maxbeta
                            beta_a(j,x) = p.mod_adapt.maxbeta;
                        elseif beta_a(j,x) < p.mod_adapt.minbeta
                            beta_a(j,x) = p.mod_adapt.minbeta;
                        end
                    end
                else
                    beta_a(:,x) = .1;
                end

                %% --- Plot Routine
%                 if p.plotinloop==1
%                     phiplot = reshape(phi_a(:,:,x),[p.simulation.ndata,p.num_terms]);
%                     %             rhoplot = reshape(rho_a(:,:),size(rho,1),size(rho,3));
%                         perttitle=RunPlotRoutinep(h,a,x,its,y0,p,colors,...
%                             sum_abs_error_a(1:x),max_corr_a(1:x),...
%                             pert_index(a,:),p_t_a(:,x),phiplot,rho_a(:,:),error);
% 
%                 end
                    
                if max_corr_a(x)/sum_abs_error_a(x) == max_corr_a(x-1)/sum_abs_error_a(x-1) 
                     p_t_a(:,x:its) = zeros(size(p_t_a,1),its-x+1);
                    phi_a(:,:,x:its) = zeros(size(phi_a,1),size(phi_a,2),its-x+1);
                    error(:,x:its) = 54321*ones(size(error,1),its-x+1);
                    sum_abs_error_a(x:its)=54321*ones(1,its-x+1);
                    max_corr_a(x:its)=0;
                break;
                end
                
            else
                p_t_a(:,x:its) = zeros(size(p_t_a,1),its-x+1);
                phi_a(:,:,x:its) = zeros(size(phi_a,1),size(phi_a,2),its-x+1);
                error(:,x:its) = 54321*ones(size(error,1),its-x+1);
                sum_abs_error_a(x:its)=54321*ones(1,its-x+1);
                max_corr_a(x:its)=0;
%                 outputs_a(x:its)
                break;
            end
        end
        
        p_t(a,:,:) = p_t_a; 
        phi(a,:,:,:) = phi_a; 
%         error(a,:,:) = error_a;
        sum_abs_error(a,:) =  sum_abs_error_a; 
        max_corr(a,:) =  max_corr_a;
        outputs(a,:) =  outputs_a;
        beta(a,:,:) = beta_a;
        gamma(a,:,:) =  gamma_a;
        rho(a,:,:) =  rho_a; 
        mu(a,:,:) = mu_a;
        
        %clearvars y0 y_all y_ave
%         end
    close 
end
% choose best from parallel results
best_error = nom_error;
best_corr = nom_corr;
continueontosecondphase=0;
for a=1:num_models
    if best_corr_p(a)/best_error_p(a) > best_corr/best_error

      best_mod=best_models_p(1,a);

        best_error = best_error_p(a);
        best_corr = best_corr_p(a);
        y_best = y_best_p(:,a);
        best_rho = best_rho_p(:,a);
%                         best_rho_hist = rho_a;
        best_models = best_models_p(:,a);
        best_pertnum = best_pertnum_p(a,:);
        best_mu = best_mu_p(:,a);
        best_beta = best_beta_p(:,a);
%                         endits = x;
%         if p.plotinloop==1
%             best_perttitle = best_perttitle_p(a);
%         else
%             best_perttitle='';
%         end
        continueontosecondphase=1;
    end
end
if continueontosecondphase
%% plot winner from first round
    if p.plotinloop
        figure;
        plot(y_best,'r'); hold on;
        plot(p.yhat,'--g');
        plot(p.Y,'b'); 
        title(['Best from initial run: $\hat{y} = ' char(best_models(1).eqn_form) '$'],'interpreter','latex');
        l=legend('$\hat{y}_0$','$\hat{y}_{end}$','$y^*$');
        set(l,'interpreter','latex','Fontsize',12);
    end
    
    
    %% Gradient descent for best perturbation case, starting from best model
    % keyboard
    
    gamma2(1:p.num_terms,1:its2)=0; %reset can_mod = nom_mod
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
    p_t2 = zeros(p.num_terms,its);
    % delM2=p_t2;
    pertlist='';
    itsstart=2;
else
    disp('No better model forms found.');
end

elseif p.continuept2
    fprintf(['Continuing iterations...\n']);
    itsstart=its2;
    its2=2*its2;
end
if continueontosecondphase
for count=1:length(best_pertnum)
    pertlist=strcat(pertlist,',',char(p.intvars(best_pertnum(count))));
end
if p.plotinloop2
    h=figure;
end
pause(1.0);
for x = itsstart:its2
    pass3=1;
    
   % p.mod_adapt.algebra = 0;
    %% --- model adaptation routine
    if x == its2
        keyboard
    end
    try
        
        [p_t2(:,x),phi2(:,:,x),error2(:,x),sum_abs_error2(x),max_corr2(x),outputs,pass3] = Sensitivity_compute(best_models,beta2(:,x-1),gamma2(:,x-1),rho2(:,x-1),p_t2(:,x-1),p);
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
        
        %% --- Update best model
        if max_corr2(x)/sum_abs_error2(x) > best_corr/best_error
    %       if max_corr(end) > best_corr
            best_mod = best_models(1);

            best_error=sum_abs_error2(x);
            best_corr = max_corr2(x);
            y_best = y02(:,1,x);
%                 y_best = y_ave2(:,1,end);
            best_rho = rho2(:,x-1);

        end

        if max_corr2(x)>win_criteria
    %         best_gamma = gamhist2(:,x);
            disp(['Correlation criteria achieved (max_corr2 = ' num2str(max_corr2(x)) ')']);
            its_FTW = x;
            break;
        end
            
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
                elseif abs(dMC2(j)) > p.mod_adapt.maxDMC
                    beta2(j,x) = beta2(j,x)-p.mod_adapt.betastep;
                    %         elseif beta(j,end) > p.mod_adapt.p.init_beta(j)
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
        if p.plotinloop2
            phiplot = reshape(phi2(:,:,x),[p.simulation.ndata,p.num_terms]);
            perttitle=RunPlotRoutine(h,1,x,its2,y_ave2,p,colors,...
                sum_abs_error2(1:x)',max_corr2(1:x)',...
                pertnum2,p_t2(:,x),phiplot,rho2,error2);
        end
        

    end
end


end 
if p.plotinloop2
        figure;
        plot(p.yhat,'--g'); hold on;
         plot(p.Y,'b');
        plot(y_best,'.-r'); 
        %hold on;
        title(['Best final model: $\hat{y} =' char(best_models(1).eqn_form) '$'],'interpreter','latex');
        l=legend('$\hat{y}_0$','$y^*$','$\hat{y}_{end}$');
        set(l,'interpreter','latex','Fontsize',12);
end
if p.save==1
    
    savefile = [p.savepath '_its' num2str(its) '-' num2str(its2)];

    save(savefile); 
%     best_mod.eqn_sym
end
out.nom_mod = p.nom_mod;
out.nom_error = nom_error;
out.nom_corr = nom_corr;
out.best_mod = best_mod;
out.best_models = best_models;
out.y_best = y_best;
out.best_error = best_error;
out.best_corr = best_corr;
out.best_rho=best_rho;
if continueontosecondphase
disp(['best model: ' char(out.best_mod.eqn_form)]);
disp(['error improvement: ' num2str((out.nom_error-out.best_error)/out.nom_error*100) '%']);
disp(['correlation improvement: ' num2str((out.best_corr-out.nom_corr)/out.nom_corr*100) '%']);
end
catch ME
    disp(['error: ' ME.message]); 
    for i=1:length(ME.stack), disp(ME.stack(i)); end
    disp('execution paused in MSAMp.m. Press F5 to exit');
    keyboard
end
end