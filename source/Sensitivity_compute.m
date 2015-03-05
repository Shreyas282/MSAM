function [p_t,phi,error,sum_abs_error,max_corr,outputs,pass] = Sensitivity_compute(models,beta,gamma,rho,p_t,p)
% this function takes a model and its perturbations, simulates their
% output, compares it to the target, and calculates the update size for the
% model, applies it to the model, and retuns the updated model form. 

% INPUT
% models:
% beta:
% gamma:
% rho: 
% p:     structure of parameter settings 
%        .mc_opt: model complexity option
%        .struct_err_opt
%        .corr_trigger: trigger to switch StructError calculation type
%        .err_trigger: ""
%        .max_beta
%        .min_beta
%        .adapt_beta
%        .betastep
%        .init_beta
%        .minDMC
%        .maxDMC
%        .mustart
%        .agg
%        .Y : target output

% OUTPUT
% [models, beta, gamma, rho, p_t, p, pass]
sum_abs_error = 0;
max_corr=0;
phi = zeros(p.simulation.ndata,size(beta,1));
error = zeros(p.simulation.ndata,1);
outputs.y0 = [];
outputs.y_all= [];
outputs.y_ave = [];
outputs.dMC = [];
if p.mod_adapt.output_dydtheta
    outputs.dydtheta=[];
    outputs.dydtheta_norm=[];
end
%% Run perturbed model simulations
% [y0, y_all, y_ave, dMC, MC] = RunPertModels_sym(models,beta,p); 
if p.mod_adapt.algebra
%     try
        RPout = RunModels_alg_sym(models,beta,p); 
%         RPout = RunModels_alg_sym_nonlinP(models,beta,p);
%     catch
%         keyboard
%     end
else
%     disp('running pert models');
    RPout = RunPertModels_sym(models,beta,p); 
%     disp('exited pert models');
end
% catch
%     fprintf(['in modeladaptationroutine.m']);
%     keyboard
% end
dMC = RPout.dMC;
MC = RPout.MC;
y0 = RPout.y0;
y_all = RPout.y_all;
y_ave = RPout.y_ave;
if p.mod_adapt.mc_opt2 
    dPhiT_dG = RPout.dPhiT_dG;
end
%% simultation failure routine
            if sum(dMC == 0)>0 % means simulation failed
                % erase past iteration and make changes more conservative
        %         beta = beta.*.5; % reduce step size
%                 mu = mu.*.25; % reduce confidence
%                 p_t = p_t;
                    
%                 gamma= gamma; % get previous gamma values before crash
%                 gamhist(:,) = gamhist(:,end-1);
%                 rho = rho;
                
%                 rho = rho-gamma; % remove last gamma step
%                 %update gamma with less confidence
%                 gamma = 0.25.*p_t.*rho;
%                 rho= rho+gamma; 

                bite=0;
            else
                bite=1;
            end
%         end
        if bite
            
      %% get NLS structural error
      [p_t,phi] = StructError(y0(:,1),y0(:,2:end),p.Y,dMC,p.mod_adapt.struct_err_opt);
 
        if (any(isnan(p_t)))
            disp('NaN values in p_t');
            keyboard
        end
    % calculate error and correlation
%         error = p.Y-y_ave(:,1);
%         sum_abs_error = sum(abs(p.Y-y_ave(:,1)));
% %         max_corr= max(xcorr(p.Y,y_ave(:,1),'coeff')); 
%         tmp = corrcoef(p.Y,y_ave(:,1));
%         max_corr = (tmp(1,2));
        error = p.Y - y0(:,1);
        sum_abs_error = sum(abs(p.Y-y0(:,1)));
%         max_corr = R2(p.Y,y0(:,1));
        tmp = corrcoef(p.Y,y0(:,1));
        max_corr = tmp(1,2);
outputs.y0 = y0;
outputs.y_all= y_all;
outputs.y_ave = y_ave;
outputs.dMC = dMC;
if p.mod_adapt.MC_NLS_plot == 1
    outputs.dydtheta=RPout.dydtheta;
    outputs.MC = MC;
end
if p.mod_adapt.output_dydtheta
    outputs.dydtheta=RPout.dydtheta;
    outputs.dydtheta_norm = RPout.dydtheta;
end
        % update settings if really close to target
%         if max_corr > p.mod_adapt.corr_trigger && sum_abs_error < p.mod_adapt.err_trigger
%             p.mod_adapt.struct_err_opt = 'errcorr';
%             p.mod_adapt.maxbeta=p.mod_adapt.maxbeta*.5;
%             beta = beta*.5;
%         end
%         %% adapt mu 
%         if p.mod_adapt.adapt_mu
%             mu = GetMu(p.Y,y_ave,p,p_t,beta);            
%         end
% %         gamma = [mu(1)*p_t(1)*rho(1);
% %                  mu(2)*p_t(2)*rho(2)];
%              gamma = mu.*p_t.*rho;
%         rho = rho+gamma;
% 
% % Update models - candidate and perturbed - to current state
% % try
% %% update model
%         models = UpdateModels_sym(models,gamma,p);
% % catch
% %     keyboard
% % end
% 
%         
% %% update beta size  
% if p.mod_adapt.adapt_beta
%      for j=1:length(dMC)
%         if abs(dMC(j)) < p.mod_adapt.minDMC
%             beta(j) = beta(j)+p.mod_adapt.betastep;
%         elseif abs(dMC(j)) > p.mod_adapt.maxDMC
%             beta(j) = beta(j)-p.mod_adapt.betastep;
% %         elseif beta(j,end) > p.mod_adapt.init_beta(j)
% %             beta(j) = beta(j) - p.mod_adapt.betastep;
%         end
%         if beta(j) > p.mod_adapt.maxbeta
%             beta(j) = p.mod_adapt.maxbeta;
%         elseif beta(j) < p.mod_adapt.minbeta
%             beta(j) = p.mod_adapt.minbeta;
%         end
%      end
% end


% outputs.y0 = y0;
% outputs.y_all= y_all;
% outputs.y_ave = y_ave;

%% estimate PARSIM
if p.wt_estimate
              
              WT = [1]; % if w=1 wavelet='DerGauss'; w=2 wavelet='Sombrero'; if w=3 wavelet='Gauss'; if w=4 wavelet='Morlet';
              %output_focus = [1];
              thresh_range=[1.5, 1.8, 2];
                           
              p_w = PARSIM(error,phi,WT,thresh_range,mu,rho);
              p_t = p_w;
              
end     

%% -- Check validity
error_hat = zeros(p.simulation.ndata,1);
error_one = [];
                   for ccount = 1:length(p_t)
                   error_one(:,ccount) = p_t(ccount).*phi(:,ccount);
                   error_hat = error_hat + error_one(:,ccount);
                   end
                                      
                   tmp = corrcoef(error,error_hat);
                   how_fit = tmp(1,2);
                   if ~isnan(how_fit)
                        p_t = p_t*how_fit;
                   end
                   
%                    if how_fit < 0.90
%                        p_t = sign(p_t);
%                    end
                   
pass=1;
        else
            pass=0;
        end
end