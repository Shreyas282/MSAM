function out = MSAM(p,its,its2)
%Model Structure Adaptation Method
% William La Cava 2015

out=[];
disp('\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\');
disp(' Running MSAM version 1.0.0-wgl');
disp('///////////////////////////////');

disp(['Number of round robin iterations: ' num2str(its)]);
disp(['Number of final adaptation iterations: ' num2str(its2)]);

if ~p.continuept2    
     colors = ['r','b','g','m','c'];
    %% initial error, correlation and string distance
    num_models = length(p.intvars)^p.num_terms; % number of models to try (hill climbs)
    %pert_index = zeros(length(p.intvars),p.num_terms);
    a = repmat({1:length(p.intvars)},1,p.num_terms);
    pert_index = allcomb(a{:})';
    nom_error = sum(abs(p.Y-p.yhat(:,1)));
    tmp = corrcoef(p.Y,p.yhat(:,1));
    nom_corr = tmp(1,2);
   
    max_corr = nom_corr;
    best_corr  = nom_corr;
    %best_gamma = [0,0];
    best_rho = zeros(p.num_terms,1);
    best_error = nom_error;
    sum_abs_error=nom_error;
    best_mod = p.nom_mod; % best model from all runs
    y_best = p.yhat;
    pass(1:num_models) = 0;
    p.mod_adapt.bestpass =1;
    p.mod_adapt.valleycontrol = 1;
    win_criteria = .9995;
       
    nom_gamma = ones(p.num_terms,1);
    phi = zeros(p.simulation.ndata,p.num_terms,its,num_models);
    error = zeros(p.simulation.ndata,its);
    gamma=zeros(p.num_terms,its,num_models); %reset can_mod = nom_mod
    rho = zeros(p.num_terms,its,num_models);
    beta=[];
    p_t = zeros(p.num_terms,its,num_models);
    p_w = zeros(p.num_terms,its,num_models);
    delM = zeros(p.num_terms,its,num_models);
    mu=[];
    chosen=struct();
   % pertnum=zeros(p.num_terms,num_models);
   % models = [];
    %% run hill climb
   
%    n=5;
 continueontosecondphase=0;
for a = 1:num_models
    %reset variables for each hill climb
    disp(['Running model adaptation ' num2str(a) ' of ' num2str(num_models)]);
    can_mod = p.nom_mod;
    beta(:,:,a)=p.init_beta;

    mu(:,:,a) = p.init_mu;
    %     nom_mod.eqn_sym = GetEqnSym(nom_mod);
    %     can_mod.eqn_str = GetEqnStr_sym(can_mod);

    % choose unique set of perturbation types
%         while pass(a) ~= 1
            [chosen,can_mod] = ChoosePertType(can_mod,p,pert_index(:,a));
            models(:,a) = [can_mod;chosen];

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
        for j=1:length(t_opts)
            % if the perturbation variable is already part of the term
            tmp=regexp([char(can_mod.terms(t_opts(j)).val)],char(p.intvars(pert_index(j,a))),'once');
            if sum(tmp)>0
                % if the perturbation variable is already part of the
                % perturbed term, start rho at the exponent of that term
                rho(j,:,a)=nom_gamma(j);
            else
                % otherwise we are perturbing with a new variable, so rho
                % is the size of the perturbation

                rho(j,:,a) = beta(j,:,a);


            end
        end
        % run gradient descent for chosen perturbation set
        %         if a==1
        if p.plotinloop
            h = figure;
        end
        %         end
        for x = 2:its
            %% --- model adaptation routine
            %% ----- compute sensitivity Phi

            %       if Produce_phi==1
%                 try

                [p_t(:,x,a),phi(:,:,x,a),error(:,x,a),sum_abs_error(x,a),max_corr(x,a),outputs(x,a),pass2(x,a)] = Sensitivity_compute(models(:,a),beta(:,x-1,a),gamma(:,x-1,a),rho(:,x-1,a),p_t(:,x-1,a),p);
                if p.mod_adapt.MC_NLS_plot == 1
                    p.dydtheta=outputs(x,a).dydtheta;
                    p.MC = outputs(x,a).MC;
                end
                %            delM = p_t;
%                 catch
%                     disp('SensitivityCompute() failed');
%                     %keyboard
%                     pass2(x,a)=0;
%                 end

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
                            %         elseif beta(j,end) > p.mod_adapt.p.init_beta(j)
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
                if p.plotinloop==1
                    phiplot = reshape(phi(:,:,x,a),[p.simulation.ndata,p.num_terms]);
                    %             rhoplot = reshape(rho(:,:,a),size(rho,1),size(rho,3));
                    try
                        perttitle=RunPlotRoutine(h,a,x,its,y_ave,p,colors,...
                            sum_abs_error(1:x,a),max_corr(1:x,a),...
                            pert_index(:,a),p_t(:,x,a),phiplot,rho(:,:,a),error(:,:,a));
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
                        best_pertnum=pert_index(:,a);
                        best_mu = mu(:,x,a);
                        best_beta = beta(:,x,a);
                        endits = x;
                        if p.plotinloop==1
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
%         end
    close 
   end

if continueontosecondphase
%% plot winner from first round
	if p.plotinloop
		figure;
		plot(y_best,'r');
		plot(p.yhat,'--g');
		plot(p.Y,'b'); hold on;
		title(['Best from initial run' best_perttitle]);
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

end
if continueontosecondphase && p.plotinloop2
        figure;
        plot(p.yhat,'--g'); hold on;
        plot(y_best,'r'); 
        plot(p.Y,'b'); %hold on;
        title(['Best final model' best_perttitle]);
        l=legend('$\hat{y}_0$','$\hat{y}_{end}$','$y^*$');
        set(l,'interpreter','latex','Fontsize',12);
end
if p.save==1
    
    savefile = [p.savepath p.simulation.sim_model '_mu' ...
        num2str(p.mod_adapt.mustart(1)*100) '_its' num2str(its) '-' num2str(its2)]

    save(savefile); 
    best_mod.eqn_sym
end
out.best_mod = best_mod;
out.y_best = y_best;
out.best_error = best_error;
out.best_corr = best_corr;
out.best_rho=best_rho;

end