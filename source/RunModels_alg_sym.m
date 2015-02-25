function output = RunModels_alg_sym(model_options,beta,p)
    
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
    % warning('off','all');
    % warning;

% if p.simulation.partitioning
% %     sim_model = p.simulation.sim_model_part;
%     x1_in = p.x1_in;
%     x2_in = p.x2_in;
% else
% %     sim_model = p.simulation.sim_model;
% end
%     x1 = p.x1_in;
%     dx1 = p.dx1_in;
%     u_in = p.u_in;
%     sim_model = p.simulation.sim_model;
    n=length(model_options);
    numperts = n-1;  
    ndata=p.simulation.ndata;
    
    options=simset('SrcWorkspace','current','DstWorkspace','current');

    
    % declare parameter values from stored values
     for count=1:length([p.cons(1,:)])
        eval([p.cons{1,count} ' = ' num2str(p.cons{2,count}) ';']);
		cons(count) = eval([p.cons{1,count} ';']);
     end
    
     % substitute parameters values into each model
     for count = 1:length(model_options)
         model_eqns(count) = subs(model_options(count).eqn_sym,...
         {p.cons{1,:}},p.cons(2,:));
     end
     
     % get lumped coefficients of model equation
     constantvals = coeffs(model_eqns(2));
     
     
    % declare new variables to represent coefficients
      A = sym('a',[getNumTerms(model_eqns(2)),1]);

     % sub new lumped constants into coefficient spots
         models = model_options;


     
     if ~strcmpi(p.mod_adapt.mc_opt,'pertsize')
         % define new lumped parameters for each model
         
         for count = 1:length(models)
             
            models(count) = getTerms(models(count),'mod',p);
             left = 1:length([models(count).terms]);
            for count2 = 1:length([models(count).terms])
                                
                 models(count).terms(count2).val = sym(subs(models(count).terms(count2).val,{p.cons{1,:}},cons));
                 models(count).terms(count2).val = simple(models(count).terms(count2).val);
                 try
                 if count==1
                     
                     constantval(count2) = coeffs(models(count).terms(count2).val);
%                      if constantval(count2)~=1
                        constants(count2) = eval(['sym(''a' num2str(count2) ''');']);
%                      end
                     
%                  end              
%                  models(count).terms(count2).val = subs(models(count).terms(count2).val,constantval,constants(count2));
                    if constantval(count2)~=1
                        models(count).terms(count2).val = subs(models(count).terms(count2).val,constantval(count2),constants(count2));
                    else
                        models(count).terms(count2).val = models(count).terms(count2).val*constants(count2);
                    end
%                  lumpedparams(count) = eval([char(constantval(count2)) ';']);
                 else
                    tmp=1;
                    while(constantval(left(tmp)) ~= coeffs(models(count).terms(count2).val))
                        tmp=tmp+1;
                    end             
                    if constantval(left(tmp))~=1
                        models(count).terms(count2).val = subs(models(count).terms(count2).val,constantval(left(tmp)),constants(left(tmp)));
                    else
                        models(count).terms(count2).val = models(count).terms(count2).val*constants(left(tmp));
                    end
                    left(tmp)=[];
                 end
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
     
    % declare variable values
    for count = 1:length(p.allvars)
        eval([char(p.allvars(count)) ' = p.' char(p.allvars(count)) ';']);
    end
    
    % define perturbation from beta
    for count = 1:length(beta)
        eval(['b' num2str(count) ' = beta(count);']);
    end
    
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
%         model_options(j).eqn_sym = subs(model_options(j).eqn_sym,{'b1','b2'},[b1 b2]);
%         model_eqns(j) = subs(model_eqns(j),{'b1','b2'},[b1 b2]);
%         outstr = GetEqnStr_sym(model_eqns(j),p.allvars,'eqn');
        model_str = regexprep(char(model_eqns(j)),'\^','\.\^');
        model_str = regexprep(model_str,'\*','\.\*');
        model_str = regexprep(model_str,'\/','\.\/');
        model_str = regexprep(model_str,'abs\(q\)','sign\(q\)\.\*abs\(q\)');
        model_str = regexprep(model_str,'abs\(qdot\)','sign\(qdot\)\.\*abs\(qdot\)');
        model_str = regexprep(model_str,'abs\(qddot\)','sign\(qddot\)\.\*abs\(qddot\)');
        model_str = regexprep(model_str,'abs\(dx1\)','sign\(dx1\)\.\*abs\(dx1\)');
        model_str = regexprep(model_str,'abs\(dx2\)','sign\(dx2\)\.\*abs\(dx2\)');
        model_str = regexprep(model_str,'abs\(dx3\)','sign\(dx3\)\.\*abs\(dx3\)');
        model_str = regexprep(model_str,'abs\(dx4\)','sign\(dx4\)\.\*abs\(dx4\)');
        model_str = regexprep(model_str,'abs\(x1\)','sign\(x1\)\.\*abs\(x1\)');
        model_str = regexprep(model_str,'abs\(x2\)','sign\(x2\)\.\*abs\(x2\)');
        model_str = regexprep(model_str,'abs\(x3\)','sign\(x3\)\.\*abs\(x3\)');
        model_str = regexprep(model_str,'abs\(x4\)','sign\(x4\)\.\*abs\(x4\)');
        try
            
            y0(:,j) = (eval([model_str ';']));
            pass=1;
        catch

            fprintf(['Model form ' model_str ...
                ' failed to output a valid response in simulink.\n']);
            keyboard
            pass=0;
        end
        y0(isnan(y0(:,j)),j)=0;
        y0(isinf(y0(:,j)),j)=0;
        error = sum(abs(p.Y-y0(:,j)));
        if isnan(error)
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
                            
                    sum_error(j)=0;
                    try
                        sim(p.simulation.sim_model,[],options);
                    catch
                        fprintf(['Model form ' char(model_options(j).eqn_sym) ...
                            ' failed to output a valid response in simulink.\n']);
%                         pause(.01);
%                         keyboard
                        %                     pass=0;
                        e(:,j,q)=(10^35)*ones(ndata,1);
                        break
                    end
					y_tmp(isnan(y_tmp))=0;
                    y_tmp(isinf(y_tmp))=0;
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
                        try
                            eval([char(A(count)) ' = dTheta(q,count);']);
                        catch
                            disp('warning: different number of terms than parameters');
                        end
                    end
                    
                    try
                         y_tmp = (eval([model_str ';']));
                    catch
                        fprintf(['Model form ' char(model_options(j).eqn_sym) ...
                            ' failed to output a valid response in simulink.\n']);
%                         pause(.01);
                        keyboard
                        %                     pass=0;
                        e(:,j,q)=(10^35)*ones(ndata,1);
                        break
                    end
                     y_tmp(isnan(y_tmp))=0;
                    y_tmp(isinf(y_tmp))=0;
                    y_all(:,j,q) = y_tmp;
                    if mod(q,2)==0 % if q is even
                        dydtheta(:,j,counter) = (y_all(:,j,q) - y_all(:,j,q-1))/(2*pert);
%                         dydtheta(:,j,counter) = (y_all(:,j,q) - y_all(:,j,q-1))/max(abs(dTheta(q,:)-dTheta(q-1,:)));
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
%                         dMC(j-1) = (sum(dydtheta_norm(:,j))/MC(1));
                        
                    end

                end
                if j>1
                    if p.mod_adapt.mc_opt2
                            dMC(j-1) = beta(j-1);
                    else
                        dMC(j-1) = (sum(dydtheta_norm(:,j))/MC(1));
                    end
                end
                MC(j) = sum(dydtheta_norm(:,j));
                

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
    if p.mod_adapt.output_dydtheta && strcmpi(p.mod_adapt.mc_opt,'pertsize')~=1
        output.dydtheta = dydtheta;
        output.dydtheta_norm = dydtheta_norm;
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