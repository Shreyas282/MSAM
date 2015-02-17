function output = RunPertModels_sym(model_options,beta,p)
    
    % INPUTS
    % true_model - target model
    % model_options - candidate model and its perturbed states
    % beta - perturbance exponent
    % dm_opt - determines how model complexity is calculated
    %        - 'var' : calculated using variance of y_all 
    %        - 'varcorr': calculated using variance of correlation of
    %        y_all with target model output
    %        - 'pardif' : calculated from the square sum of the partial
    %        differentials of y with respect to m, c, and k
    
    % OUTPUTS
    % output.
    % Y - target output
    % y0 - model options sim output with true parameters
    % y_all - all model options outputs with all parameter values
    % y_ave - model options outputs averaged over parameter values
    % MC - model complexity
    % dydtheta - parameter sensitivities
%     warning('off','all');
%     warning;

% if p.simulation.partitioning
% %     sim_model = p.simulation.sim_model_part;
    % x1_in = p.x1_in;
    % x2_in = p.x2_in;
% else
% %     sim_model = p.simulation.sim_model;
% end
    sim_model = p.simulation.sim_model;
    n=length(model_options);
    numperts = n-1;  
    input_set = p.simulation.input_set;% 1:impulse; 2:step; 3:random
    time_span=p.simulation.time_span;
    samp_time=time_span/(p.simulation.ndata-1);
    ndata=p.simulation.ndata;
    IC = p.simulation.IC;
    noise_on = 0;
    
    options=simset('SrcWorkspace','current','DstWorkspace','current');

    
    % declare and substitute parameter values from stored values
     for count=1:length([p.cons(1,:)])
        eval([p.cons{1,count} ' = ' num2str(p.cons{2,count}) ';']);
        cons(count) = eval([p.cons{1,count} ';']);
     end
    
%      for count = 1:length(model_options)
%              model_eqns(count) = subs(model_options(count).eqn_sym,{p.cons{1,:}},cons);
%      end
     models = model_options;
     
     if ~strcmpi(p.mod_adapt.mc_opt,'pertsize')
         % define new lumped parameters for each model
         
         for count = 1:length(models)
            models(count) = getTerms(models(count),'mod',p);
             
            for count2 = 1:length([models(count).terms])
                                
                 models(count).terms(count2).val = sym(subs(models(count).terms(count2).val,{p.cons{1,:}},cons));
                 models(count).terms(count2).val = simple(models(count).terms(count2).val);
                 try
                 if count==1
                     
                     constantval(count2) = coeffs(models(count).terms(count2).val);
                     if constantval(count2)~=1
                        constants(count2) = eval(['sym(''a' num2str(count2) ''');']);
                     end
                     
                 end              
%                  models(count).terms(count2).val = subs(models(count).terms(count2).val,constantval,constants(count2));
                    if constantval(count2)~=1
                        models(count).terms(count2).val = subs(models(count).terms(count2).val,constantval,constants);
                    end
%                  lumpedparams(count) = eval([char(constantval(count2)) ';']);
                catch
%                          keyboard
                end
                
            end
            models(count).eqn_sym = GetEqnSym(models(count));
            models(count).eqn_str = GetEqnStr_sym(models(count),p.allvars);
         end
         
            model_eqns = [models.eqn_sym];
            removals=constants==0;
            constants(removals)=[];
            constantval(removals)=[];
            lumpedparams=double(constantval);
            for cnt=1:length(constantval)
                    try
                    eval([char(constants(cnt)) ' = ' num2str(constantval(cnt)) ';']);
                    catch
                        try 
                            eval([char(constants(cnt)) ' = ' char(constantval(cnt)) ';']);
                        end
                    end
            end
%             keyboard
%          for count=1:length(constantvals)
%              lumpedparams(count) = eval([char(constantvals(count)) ';']);
%              eval([char(constants(count)) ' = ' num2str(lumpedparams(count)) ';']);
%          end
         
%          constantvals = coeffs(model_eqns(1));
% 
%          for count2 = 1:getNumTerms(model_eqns(1))
%              constants(count2) = eval(['sym(''a' num2str(count2) ''');']);
%          end
% 
%          for count = 1:length(model_options)
%             model_eqns(count) = subs(model_eqns(count),constantvals,constants);
%          end
% 
%          for count=1:length(constantvals)
%              lumpedparams(count) = eval([char(constantvals(count)) ';']);
%              eval([char(constants(count)) ' = ' num2str(lumpedparams(count)) ';']);
%          end
     else
         model_eqns = [models.eqn_sym];
     end
     
    % define perturbation from beta
    for count = 1:length(beta)
        eval(['b' num2str(count) ' = beta(count);']);
    end
%     b1 = beta(1); b2 = beta(2);
%     outstr = true_model{1};
%     sim(p.simulation.sim_model,[],options);
%     Y = y_tmp;
%     input(:,1) = u(1,1,:);
    
   MC=0;
    dMC=0;
    sum_error(1:n) = 0;
    
    % set up parameters
    if strcmpi(p.mod_adapt.mc_opt,'var') || strcmpi(p.mod_adapt.mc_opt,'varcorr')
        for count=1:length([p.cons(1,:)])
            dTheta(:,count) = [eval([p.cons{1,count}])*[.8:.2:1.2]'];
        end
        dTheta=[shuffle(dTheta(:,1)),shuffle(dTheta(:,2)),shuffle(dTheta(:,3))];
    elseif strcmpi(p.mod_adapt.mc_opt,'pardif')
        pert = .1;
        dTheta = repmat(lumpedparams,[length(lumpedparams)*2,1]);
        inc=1;
        for count=1:2:size(dTheta,1)
        % perturb each model parameter larger and smaller    
            dTheta(count:count+1,inc) = [(1+pert) 0; 0 (1-pert)]*dTheta(count:count+1,inc);
            inc=inc+1;
        end
    else%standard NLS approach
        pert=beta;
        dTheta=0;
                  
    end
    y_all(ndata,length(model_options),size(dTheta,1)) = 0;
    
    
    for j = 1:length(model_options)
        for bs=1:length(beta)
        model_options(j).eqn_sym = subs(model_options(j).eqn_sym,{['b' num2str(bs)]},eval(['b' num2str(bs)]));
        model_eqns(j) = subs(model_eqns(j),{['b' num2str(bs)]},eval(['b' num2str(bs)]));
        end
        
        outstr = GetEqnStr_sym(model_eqns(j),p.allvars,'eqn');
        
        
        try
            sim(p.simulation.sim_model,[],options);
            y0(:,j)=y_tmp;
            pass=1;
            % sim(p.simulation.sim_model);
        catch ME
            ME
            %         y0=ones(128,1)*rand*10000;
            fprintf(['Model form ' outstr ...
                ' failed to output a valid response in simulink.\n']);
%             keyboard
            pass=0;
        end
        if pass==1
      %    y0(:,j)=y_tmp;  
            e = [];
            if strcmpi(p.mod_adapt.mc_opt,'var') || strcmpi(p.mod_adapt.mc_opt,'varcorr')
                for q = 1:size(dTheta,1) 

                    for count=1:length(dTheta(1,:))
                        eval([p.cons{1,count} ' = dTheta(q,count);']);
                    end
                            
%                     m=dTheta(q,1);
%                     c=dTheta(q,2);
%                     k=dTheta(q,3);

                    sum_error(j)=0;
                    try
                        sim(p.simulation.sim_model,[],options);
                    catch
                        fprintf(['Model form ' outstr ...
                            ' failed to output a valid response in simulink.\n']);
%                         pause(.01);
%                         keyboard
                        %                     pass=0;
                        e(:,j,q)=(10^35)*ones(ndata,1);
                        break
                    end
                    y_all(:,j,q) = y_tmp;
%                     sum_y = sum(y_tmp);

%                     e(:,j,q) = Y-y_tmp;
                    if strcmpi(p.mod_adaptmc_opt,'varcorr')
                        y_corr(:,j,q) = xcorr(y_all(:,1,q),y_tmp);
                    end
                    %
%                     sum_error(j) = sum_error(j)+sum(abs(e(:,j,q)));

                end
            elseif strcmpi(p.mod_adapt.mc_opt,'pardif')
                % calculated square sum of partial differentials of y with respect to m, c, and k
                % calculate y for m*.9,m*1.1, c*.9,c*1.1, k*.9,k*1.1; use
                % y0 value for center value
                % calculate dydm, dydc, dydk 
                % sum their squares
                
                      counter=1;
                for q = 1:length(dTheta(:,1))
                    for count=1:length(dTheta(1,:))
                        eval([char(constants(count)) ' = dTheta(q,count);']);
                    end
                    
                    try
                        sim(p.simulation.sim_model,[],options);
                    catch
                        fprintf(['Model form ' outstr ...
                            ' failed to output a valid response in simulink.\n']);
%                         pause(.01);
%                         keyboard
                        %                     pass=0;
                        e(:,j,q)=(10^35)*ones(ndata,1);
                        break
                    end
                    y_all(:,j,q) = y_tmp;
                    if mod(q,2)==0 % if q is even
                        dydtheta(:,j,counter) = (y_all(:,j,q) - y_all(:,j,q-1))/(2*pert);
                        counter=counter+1;
                    end
                        
                end
               
                
                for count=1:ndata
                    if j==1
                        dydtheta_norm(count,j) = norm(reshape(dydtheta(count,j,:),[1 size(dydtheta,3)]));
                    else
                        dydtheta_norm(count,j) = norm(...
                        reshape(dydtheta(count,j,:),[1 size(dydtheta,3)])-...
                        reshape(dydtheta(count,1,:),[1 size(dydtheta,3)]));
                        dMC(j-1) = (sum(dydtheta_norm(:,j))/MC(1));
                    end

                end
                
                MC(j) = sum(dydtheta_norm(:,j));
                
%                     dydm(:,j) = (y_all(:,j,2) - y_all(:,j,1))/(2*pert);                    
%                     dydc(:,j) = (y_all(:,j,4) - y_all(:,j,3))/(2*pert);                    
%                     dydk(:,j) = (y_all(:,j,6) - y_all(:,j,5))/(2*pert);
                    
%                     dydtheta(j) = sum(dydm(:,j).^2+dydc(:,j).^2+dydk(:,j).^2);
%                     dydtheta(j) = sqrt(sum(dydm(:,j).^2+dydc(:,j).^2+dydk(:,j).^2));
%                     MC(j) = sum(sqrt(dydm(:,j).^2+dydc(:,j).^2+dydk(:,j).^2));
%                     dydtheta(j) = sum(sqrt(dydm(:,j).^2+dydc(:,j).^2+dydk(:,j).^2));
%                     dydtheta(j) = sum((abs(dydm(:,j))+abs(dydc(:,j))+abs(dydk(:,j))));
                    
%                     MC(j) = dydtheta(j);
                   
%                 if j>1
% %                         dMC(j-1) = (MC(j)-MC(1))/MC(1);
% %                         dMC(j-1) = abs(MC(j)-MC(1))/MC(1);
% %                     dMC(j-1) = sqrt(sum([(dydm(:,j)-dydm(:,1)).^2 + (dydc(:,j)-dydc(:,1)).^2 + (dydk(:,j)-dydk(:,1)).^2]))/dydtheta(1);
%                     for count=1:ndata
%                         dydtheta_norm(count,j) = norm(...
%                         reshape(dydtheta(count,j,1:end),[1 size(dydtheta,3)])-...
%                         reshape(dydtheta(count,1,1:end),[1 size(dydtheta,3)]));
%                     end
%                     dMC(j-1) = sum(dydtheta_norm(:,j))/MC;
% %                     dMC(j-1) = sum(sqrt([(dydm(:,j)-dydm(:,1)).^2 + (dydc(:,j)-dydc(:,1)).^2 + (dydk(:,j)-dydk(:,1)).^2]))/MC(1);
% %                     dMC(k-1) = sum([abs(dydm(:,k)-dydm(:,1)) + abs(dydc(:,k)-dydc(:,1)) + abs(dydk(:,k)-dydk(:,1))]);
% %                     if dMC(k-1) == 0 
% %                          keyboard
% %                     end
%                 end
%                     figure;
%                     plot(reshape(y_all(:,j,:),[ndata 6]));
%                     legend('m*(1-pert)',...
%                           'm*(1+pert)',...
%                           'c*(1-pert)',...
%                           'c*(1+pert)',...
%                           'k*(1-pert)',...
%                           'k*(1+pert)');
%                       title(['model ' num2str(j)]);
            elseif strcmpi(p.mod_adapt.mc_opt,'pertsize') % standard NLS
                y_all(:,j) = y0(:,j);
                if j>1
                dMC(j-1) = pert(j-1);
                end
            end
            
        
        else
            y0=0; y_all=0; y_ave=0; dMC=0; MC=0;
            output.y0 = y0;
            output.y_all= y_all;
            output.y_ave = y_ave;
            output.dMC = dMC;
            output.MC = MC;
            if p.mod_adapt.output_dydtheta
                output.dydtheta = 0;
            end
            return;
        end                  
    
        
            if strcmpi(p.mod_adapt.mc_opt,'var')
                total_scatter_time(1:n)=0;
                Scatter_t = zeros(ndata,n);
            end
            for t=1:ndata
                 y_ave(t,j) = mean(y_all(t,j,:));
                 if strcmpi(p.mod_adapt.mc_opt,'var')
                    Scatter_t(t,j) = var(y_all(t,j,:));
                    total_scatter_time(j) = total_scatter_time(j) + Scatter_t(t,j);
                 end
            end  
            if strcmpi(p.mod_adapt.mc_opt,'var')
                MC(j) = total_scatter_time(j);
            elseif strcmpi(p.mod_adapt.mc_opt,'varcorr')

                total_y_corr_scat(1:n)=0;
                y_corr_scat(2*ndata-1,n)=0;
                for t=1:2*ndata-1
                    y_corr_scat(t,j) = var(y_corr(t,j,:));
                    total_y_corr_scat(j) = total_y_corr_scat(j)+y_corr_scat(t,j);
                end
                MC(j) = total_y_corr_scat(j);

            end
    end
    output.y0 = y0;
    output.y_all= y_all;
    output.y_ave = y_ave;
    output.dMC = dMC;
    output.MC = MC;
    if p.mod_adapt.output_dydtheta
        output.dydtheta = dydtheta;
    end
%             for k=2:length(model_options)
%                     dMC(k-1) = (MC(k)-MC(1))/MC(1);
%                     
% %                     dMC(k-1) = sum(sum([(dydm(:,k)-dydm(:,1)).^2 + (dydc(:,k)-dydc(:,1)).^2 + (dydk(:,k)-dydk(:,1)).^2],2));
% %                     dMC(k-1) = sum([abs(dydm(:,k)-dydm(:,1)) + abs(dydc(:,k)-dydc(:,1)) + abs(dydk(:,k)-dydk(:,1))]);
% %                     if dMC(k-1) == 0 
% %                          keyboard
% %                     end
% %                     dMC(k-1) = 
%             end
    
%     if ~strcmpi(p.mc_opt,'varcorr')
%         for k=2:length(model_options)
%             dMC(k-1) = MC(k) - MC(1);
%         end
%     else
% %         [t,cands,iters]=size(y_all);
%             for a=1:iters                
%                 corr_inc(a,1,r) = max(xcorr(y_all(:,2,a),y_all(:,1,a)));
%                 corr_inc(a,2,r) = max(xcorr(y_all(:,3,a),y_all(:,1,a)));
%                 corr_inc(a,3,r) = max(xcorr(y_all(:,4,a),y_all(:,1,a)));
%             end
%             dMC(1,r) = var(corr_inc(:,1,r))/var(corr_inc(:,3,r));
%             dMC(2,r) = var(corr_inc(:,2,r))/var(corr_inc(:,3,r));
    
%     [t,cands,iters]=size(y_all);
%     for a=1:iters
%         corr_inc(a,1) = corr(y_all(:,2,a),y_all(:,1,a));
%         corr_inc(a,2) = corr(y_all(:,3,a),y_all(:,1,a));
%      end
%     dM_corr(inc,1) = var(corr_inc(:,inc,1));
%     dM_corr(inc,2) = var(corr_inc(:,inc,2));
%     pause
% keyboard
end