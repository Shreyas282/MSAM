% This is an umbrella program for IS-LM model components
clear all
clear

Save=0;
load('Abhijit.mat')

begin = 1; 
begin_year=1958+begin;
finish = 30;
end_year = 1958+finish;

%Model = 'consumer spending'; % choices are 'consumer spending',
%'investment', 'IS dynamics', 'money demand', 'LM_dynamics'
Model = 'IS dynamics';

Simulation = 1;
% if strcmpi(Model,'IS dynamics') || strcmpi(Model,'LM dynamics')
%     Simulation = 1;
% end

if strcmpi(Model,'consumer spending')
    savefile = ['Results/Consumer_spending_' num2str(begin_year) '_' num2str(end_year) '.mat']; 
    % C(y-tax(y)); inputs are y: GDP and taxes:
    
    u(begin:finish,1) = A(begin:finish,6); % Column G: GDP
    u(begin:finish,2) = A(begin:finish,11); % Column L: Tax Revenue
    y = A(begin:finish,7); % Column H: Personal Consumption
    
    yhat = 'a1 + a2*u(:,1) + a3*u(:,2)';
    model_str = regexprep(yhat,'\^','\.\^');
    model_str = regexprep(model_str,'\*','\.\*');
    model_str = regexprep(model_str,'\/','\.\/');
    
elseif strcmpi(Model,'investment')
    savefile = ['Results/Investment_' num2str(begin_year) '_' num2str(end_year) '.mat']; 
    % I(y,r); inputs are y: GDP, r: nominal interest rate, and tau: tax
    % revenue; the output is the net domestic investment
    
    u(begin:finish,1) = A(begin:finish,6); % Column G: GDP
    u(begin:finish,2) = A(begin:finish,19); % Column T: Nominal Interest Rate
    u(begin:finish,3) = A(begin:finish,11); % Column L: Tax Revenue
    y = A(begin:finish,17); % Column R: Net Domestic Investment
    
    yhat = 'a1 + a2*u(:,1) + a3*u(:,2)+a4*u(:,3)';
    model_str = regexprep(yhat,'\^','\.\^');
    model_str = regexprep(model_str,'\*','\.\*');
    model_str = regexprep(model_str,'\/','\.\/');
    
elseif strcmpi(Model,'IS dynamics') && Simulation==1
    savefile = ['Results/IS_dynamics_sim_' num2str(begin_year) '_' num2str(end_year) '.mat']; 
    % ydot = alpha(e-y); where e: is the government expenditure and y is
    % the GDP
    y = A(begin:finish,6); % Column G: GDP
    tmax = (length(y));
    dt = 1; % 1 year between data points
    t = 1:dt:tmax;
    %     in1.time = t;
    %     in1.signals.dimentions = 1;
    %  %   in1(begin:finish,2) = A(begin:finish,10); % Column K: Government Expenditure
    %     in1.signals.values = A(begin:finish,10); % Column K: Government Expenditure
    u(:,1) = t; 
    u(:,2) = A(begin:finish,10); % Column K: Government Expenditure
    alpha1 = 0.1;
    t1 = t';
    Init(:,1) = t1;
    Init(:,2) = repmat(y(1),size(t1));
    ydot = 'a1*(u(1) - u(2))'; % u(1): E, u(2): y
%     ydot = 'u(1)*a1*abs(u(1))^(15315970234811473/288230376151711744) - (u(2)*a1)/abs(u(1))^(378661661794032277/288230376151711744)';
%     ydot = 'u(1)*a1*abs(u(1))^(1928812644929541/36028797018963968) - (u(2)*a1)/abs(u(1))^(174759317463501585/144115188075855872)';
    model_str = regexprep(ydot,'\^','\.\^');
    model_str = regexprep(model_str,'\*','\.\*');
    model_str = regexprep(model_str,'\/','\.\/');
    
elseif strcmpi(Model,'IS dynamics') && Simulation==0
    savefile = ['Results/IS_dynamics_alg_' num2str(begin_year) '_' num2str(end_year) '.mat']; 
    % ydot = alpha(e-y); where e: is the government expenditure and y is
    % the GDP
    u(:,1) = A(begin:finish,10); % Column K: Government Expenditure
    y = A(begin:finish,6); % Column G: GDP
    tmax = (length(y));
    dt = 1; % 1 year between data points
    t = 1:dt:tmax;
    t1 = t';
    Init(:,1) = t1;
    Init(:,2) = repmat(y(1),size(t1));
    [YO,goodness] = fit(t1,y,'cubicinterp');
   %[YO,goodness] = fit(t1,y,'smoothingspline');
   %[YO,goodness] = fit(t1,y,'poly2');
    [ydot,yddot] = differentiate(YO,t1);
%     yhat = 'u(:,1) - 1/a1*ydot(:,1)+ Init(:,2)'; % Init(:,2): y0, 
    yhat = 'Init(:,2) + u(:,1)*abs(ydot(:,1))^(0.1738) - ydot(:,1)/(a1*abs(ydot(:,1))^(0.2166))';
    model_str = regexprep(yhat,'\^','\.\^');
    model_str = regexprep(model_str,'\*','\.\*');
    model_str = regexprep(model_str,'\/','\.\/');
    
elseif strcmpi(Model,'money demand')
    savefile = ['Results/Money_demand_' num2str(begin_year) '_' num2str(end_year) '.mat']; 
    % m_d = L(y,r,i); inputs are y: GDP, r: interest rate, and i: inflation
    % rate
    
    u(begin:finish,1) = A(begin:finish,6); % Column G: GDP
    u(begin:finish,2) = A(begin:finish,19); % Column T: Nominal Interest Rate
    u(begin:finish,3) = A(begin:finish,20); % Column L: inflation rate
    y = A(begin:finish,4); % Column E: M1
    
    yhat = 'a1 + a2*u(:,1) + a3*u(:,2) + a4*u(:,3)';
    model_str = regexprep(yhat,'\^','\.\^');
    model_str = regexprep(model_str,'\*','\.\*');
    model_str = regexprep(model_str,'\/','\.\/');
    
elseif strcmpi(Model,'LM dynamics') && Simulation==1
    savefile = ['Results/LM_dynamics_sim_' num2str(begin_year) '_' num2str(end_year) '.mat']; 
    % rdot = beta(L(y,r,i)-i)
    
    y = A(begin:finish,19); % Column T: Nominal interest rate
    tmax = (length(y));
    dt = 1; % 1 year between data points
    t = 1:dt:tmax;
    %     in1.time = t;
    %     in1.signals.dimentions = 1;
    %  %   in1(begin:finish,2) = A(begin:finish,10); % Column K: Government Expenditure
    %     in1.signals.values = A(begin:finish,10); % Column K: Government Expenditure
    u1(:,1) = t;
    u1(:,2) = A(begin:finish,4); % Column E: M1 
    u2(:,1) = t;
    u2(:,2) = A(begin:finish,20); % Column U: Inflation Rate
    a1 = 0.2;
    a2 = 0.2;
    t1 = t';
    Init(:,1) = t1;
    Init(:,2) = repmat(y(1),size(t1));
    ydot = 'a1*u(1) - a2*u(2) - a3*u(3)';
    model_str = regexprep(ydot,'\^','\.\^');
    model_str = regexprep(model_str,'\*','\.\*');
    model_str = regexprep(model_str,'\/','\.\/');
    
elseif strcmpi(Model,'LM dynamics') && Simulation==0
    savefile = ['Results/LM_dynamics_alg_' num2str(begin_year) '_' num2str(end_year) '.mat']; 
    % ydot = alpha(e-y); where e: is the government expenditure and y is
    % the GDP
    u(:,1) = A(begin:finish,4); % Column E: M1 
    u(:,2) = A(begin:finish,20); % Column U: Inflation Rate
    y = A(begin:finish,19); % Column T: Nominal interest rate
    tmax = (length(y));
    dt = 1; % 1 year between data points
    t = 1:dt:tmax;
    t1 = t';
    Init(:,1) = t1;
    Init(:,2) = repmat(y(1),size(t1));
   % [YO,goodness] = fit(t1,y,'cubicinterp');
   [YO,goodness] = fit(t1,y,'smoothingspline');
   %[YO,goodness] = fit(t1,y,'poly2');
    [ydot,yddot] = differentiate(YO,t1);
    yhat = 'a1/a3*u(:,1) - a2/a3*u(:,2) - 1/a3*ydot(:,1) + Init(:,2)';
    model_str = regexprep(yhat,'\^','\.\^');
    model_str = regexprep(model_str,'\*','\.\*');
    model_str = regexprep(model_str,'\/','\.\/');
end


mu=0.1;
%mu=1;
max_iteration=20;
output_focus=[1];
%cascade=1;
make_movie=0;
enableBreak=1;
Continuous_plot=0;

if strcmpi(Model,'consumer spending') 
    par_start = [1 1 1];
    par_num = 1:3;
elseif  strcmpi(Model,'investment')
    par_start = [1 1 1 1];
    par_num = 1:4;
elseif strcmpi(Model,'IS dynamics') && Simulation==1
    par_start = [0.2];
    par_num = [1];
elseif strcmpi(Model,'IS dynamics') && Simulation==0
    par_start = [0.1];
    par_num = [1];
elseif  strcmpi(Model,'money demand')
    par_start = [1 1 1 1];
    par_num = 1:4;
elseif strcmpi(Model,'LM dynamics') && Simulation==1
    par_start = [0.01 0.1 0.1];
    par_num = [1:3];
elseif strcmpi(Model,'LM dynamics') && Simulation==0
    par_start = [0.01 1 1];
    par_num = [1:3];
end
par = par_start ;
j=length(par_start);

i=1;

for ii=1:max_iteration
    
    %parameters(:,incremental_adapt)=par_in';
    pp(1:length(par_start),i) = par';   
%     DTheta(:,i)=(par_true(1:length(par))'-pp(1:length(par),i))';
    
if strcmpi(Model,'consumer spending') 
    p.cons = {'a1', 'a2', 'a3'; 
           par(1),  par(2), par(3)};
    for z=1:size(p.cons,2)
        eval([p.cons{1,z} '= p.cons{2,z};']);
    end
    y_tmp = (eval([model_str ';']));
    
elseif strcmpi(Model,'investment')
    p.cons = {'a1', 'a2', 'a3' 'a4'; 
           par(1),  par(2), par(3), par(4)};
    for z=1:size(p.cons,2)
        eval([p.cons{1,z} '= p.cons{2,z};']);
    end
    y_tmp = (eval([model_str ';']));
elseif strcmpi(Model,'IS dynamics') && Simulation==1
    p.cons = {'a1'; 
    par(1)};
    for z=1:size(p.cons,2)
        eval([p.cons{1,z} '= p.cons{2,z};']);
    end
    sim('IS_dyn');
    
elseif strcmpi(Model,'IS dynamics') && Simulation==0    
    p.cons = {'a1'; 
    par(1)};
    for z=1:size(p.cons,2)
        eval([p.cons{1,z} '= p.cons{2,z};']);
    end
    y_tmp = (eval([model_str ';']));
    
 elseif strcmpi(Model,'money demand')
    p.cons = {'a1', 'a2', 'a3' 'a4'; 
           par(1),  par(2), par(3), par(4)};
    for z=1:size(p.cons,2)
        eval([p.cons{1,z} '= p.cons{2,z};']);
    end
    y_tmp = (eval([model_str ';'])); 
    
    elseif strcmpi(Model,'LM dynamics') && Simulation==1
    p.cons = {'a1' 'a2' 'a3'; 
    par(1) par(2) par(3)};
    for z=1:size(p.cons,2)
        eval([p.cons{1,z} '= p.cons{2,z};']);
    end
    sim('LM_dyn');
    
elseif strcmpi(Model,'LM dynamics') && Simulation==0    
    p.cons = {'a1' 'a2' 'a3'; 
    par(1) par(2) par(3)};
    for z=1:size(p.cons,2)
        eval([p.cons{1,z} '= p.cons{2,z};']);
    end
    y_tmp = (eval([model_str ';']));
end
    
%     A = par(1);
%     eps = par(2);

    inc=1;
    for kk = output_focus;
        y1(:,inc) = y_tmp(:,kk);
        inc=inc+1;
    end
    
    e0(:,:,i)=y-y1; % this is the simulation error
    mean_error(i,:) = mean(abs(y-y1));
%     e0Plot = 1;
%     
%     if (e0Plot == 1)
%         figure(2);
%                hold off;
%         plot(e0(:,1,i)); hold on;
%         if (ii>1)
%         plot(Phi_t*DTheta(par_num,i),'r')
%         leg=legend('$\epsilon_0$','$\Phi_t \Delta \Theta$');
%         set(leg,'interpreter','latex');
%         else
%             max_error = max(abs(e0(:,1,i)));
%             leg=legend('$\epsilon_0$');
%             set(leg,'interpreter','latex');
%         end
% %         ylim([0 max_error]);
% %         hold on
% %         subplot(2,1,2), plot(e0(:,2,i))
% %         hold on
%         drawnow
%         
%     end
    
    pp(length(par)+1,i)=sum(sum(abs(e0(:,:,i))));
    par_error=0;
%     for ip=par_num
%         par_error =par_error + ((par_true(ip) - par(ip))/par(ip))^2;
%     end
    pp(length(par)+2,i)=par_error;
    
    pert = 0.01;
    
    delta=pert*eye(length(par));
    for l=1:j
        dpar(:,l) = par' + delta(:,par_num(l))*par(par_num(l));
    end
    
    for l=1:j;
        
%         p.cons = {'mu1', 'mu2', 'mu3', 'g1', 'g2'; 
%            dpar(1,l),  dpar(2,l),  dpar(3,l), dpar(4,l), dpar(5,l)};
        
        if strcmpi(Model,'consumer spending') 
            a1 = dpar(1,l);
            a2 = dpar(2,l);
            a3 = dpar(3,l);
            y_tmp = (eval([model_str ';']));
        elseif strcmpi(Model,'investment')
            a1 = dpar(1,l);
            a2 = dpar(2,l);
            a3 = dpar(3,l);
            a4 = dpar(4,l);
            y_tmp = (eval([model_str ';']));
        elseif strcmpi(Model,'IS dynamics') && Simulation==1
            a1=dpar(:,l);
         %   a2=dpar(:,l);
            sim('IS_dyn');
        elseif strcmpi(Model,'IS dynamics') && Simulation==0
            a1 = dpar(l,1);
            y_tmp = (eval([model_str ';']));
        elseif strcmpi(Model,'investment')
            a1 = dpar(1,l);
            a2 = dpar(2,l);
            a3 = dpar(3,l);
            a4 = dpar(4,l);
            y_tmp = (eval([model_str ';']));
         elseif strcmpi(Model,'money demand')
            a1 = dpar(1,l);
            a2 = dpar(2,l);
            a3 = dpar(3,l);
            a4 = dpar(4,l);
            y_tmp = (eval([model_str ';']));
          elseif strcmpi(Model,'LM dynamics') && Simulation==1
            a1 = dpar(1,l);
            a2 = dpar(2,l);
            a3 = dpar(3,l);
            sim('LM_dyn');
        elseif strcmpi(Model,'LM dynamics') && Simulation==0
            a1 = dpar(l,1);
            a2 = dpar(2,l);
            a3 = dpar(3,l);
            y_tmp = (eval([model_str ';']));  
        end
        inc=1;
        y2=[];
    for kk = output_focus;
        y2(:,inc) = y_tmp(:,kk);
        inc=inc+1;
    end
        dele(:,:,l) = (y2-y1)/(pert);
        if isnan(any(dele))
            keyboard
        end
    end
    
 %   a1=p.cons{2,1}; a2=p.cons{2,2}; a3=p.cons{2,3}; 
    
    dele1=dele;
    
    
%     esim(:,:,i)=zeros(size(y0,1),length(output_focus));
    
%     for k = 1:length(output_focus);
%         for l=1:j
% %             del_par(l) = (par_true(par_num(l))-par(par_num(l)))/par_true(par_num(l));
%        %     del_par(l) = (par_true(par_num(l))-par(par_num(l)));
%             esim(:,k,i) = esim(:,k,i) + del_par(l)*dele1(:,k,l);
%         end
%     end
    
      e0Plot = 1;
    
    if (e0Plot == 1)
        figure(1);
               
%        subplot(121); hold off;
      if ii==1, plot(reshape(y,size(y,1)*size(y,2),1),'b');  end
      hold on;
       plot(reshape(y1,size(y1,1)*size(y1,2),1),'--','color',[ii/max_iteration 0 0]);
%        subplot(122); hold off;
%         plot(e0(:,1,i),'r'); hold on;
        %if (ii>1)
        %plot(Phi_t*DTheta(par_num,i),'r')
%         plot(esim(:,1,i))
%         leg=legend('$y^*$','$\hat{y}$');
%         set(leg,'interpreter','latex');
%         else
%             max_error = max(abs(e0(:,1,i)));
%             leg=legend('$\epsilon_0$');
%             set(leg,'interpreter','latex');
%         end
%         ylim([0 max_error]);
%         hold on
%         subplot(2,1,2), plot(e0(:,2,i))
        hold on
        drawnow
        
    end
    
    for k=1:length(output_focus)
        dele3(:,:) = dele1(:,k,:);
        ct = corrcoef([dele3(:,:)]);
        ct_record(:,:,k,i) = ct;
    end
    
    %if cascade==1
        
        e0_vect = mat_2_vect(e0(:,:,i));
        e0_vect = e0_vect';
        for l=1:j
            dele2(:,l)=mat_2_vect(dele1(:,:,l));
        end
        outputs=1;
    %end
    
    %Gauss-Newton
    %if cascade==1
        
        Phi_t = [dele2];
        Phi_t(isnan(Phi_t)) = 0;
            Phi_t(isinf(Phi_t)) = 0;
%         try
            testmat = Phi_t'*Phi_t;
            phi_index(ii) = sqrt(min(eig(testmat))/max(eig(testmat)));
            phi_CN(ii) = max(svd(testmat))/min(svd(testmat));
%         catch
%             keyboard
            
%         end
        p_t0 = lscov(Phi_t,e0_vect);
        par_t = (p_t0.*par(par_num)');
    %end
    

    for l=1:j
        par(par_num(l))=par(par_num(l))+mu*(par_t(l));
        par_trecord(l,i) = par_t(l);
    end
    
%     sigDisplay = 3;
%     if make_movie ==1
%         if i == 1
%             blah1 = max(max(coef_mat(:, :,sigDisplay)));
%             blah2 = min(min(coef_mat(:, :,sigDisplay)));
%             
%         end
%         figure(1)
%         mesh((fliplr(sig(:,:,2).*coef_mat(:, :, sigDisplay)))')
%         axis([1 n 1 limit2 blah2 blah1])
%         parm_frame(i) = getframe;
%     end
    
    tmpPar = par';
    i = i+1;

%     save(loadfile, 'i', 'p', 'par', 'par_true', 'pp', 'e0', 'ct_record',  'par_trecord', 'mu', '-mat');

    
    par_error = 0.;
%     for ip=1:j
%         par_error =par_error+ ((par_true(par_num(ip)) - par(ip))/par(par_num(ip)))^2;
%     end
    
%    disp(['mean error: ' num2str(mean(mean_error(end,:)))]);
   
%    pp(length(par)+2,i)=par_error;
%    pp(5,i);
    
    if enableBreak==1
         if (mean(mean_error(end,:))<10^(-5)) && (enableBreak ==1)
        % if (pp(4,i-1) < 0.0001) && (enableBreak ==1)
           break;
         elseif (ii>1)
%              dif = mean_error(end)-mean_error(end-1)
                 if (abs(mean(mean_error(end,:))-mean(mean_error(end-1,:)))<10^(-8))
                     break;
                 end
        end
     end
    if Continuous_plot==1
        figure(1);
        [row,col] = size(pp);
%         par_array = zeros(row,col);
               
        
%         for kp=1:length(par_true)
%                  par_array(kp,1:length(pp)) = par_true(kp);
%         end
        
        for kp = 1:row
            subplot(row,1,kp);
            plot(pp(kp,:),'x')
            hold on
%             plot(par_array(kp,:),'r')
%             plot(repmat(par_true(kp),[1 length(pp)]),'r')
            if kp <= row-2
%                 plot(repmat(par_true(kp),[1 length(pp)]),'r')
                ylabel(['Par' num2str(kp)])
            else
                plot(repmat(pp(kp),[1 length(pp)]),'r');
                ylabel('Error');
            end
            drawnow
        end
    end
    
end % end of max_iteration


if Save==1
    if strcmpi(Model,'IS dynamics') && Simulation==1
        save(savefile, 'ydot', 'pp', 'u', 'y', 'y1', 'mu', '-mat');
    
    elseif strcmpi(Model,'LM dynamics') && Simulation==1
        save(savefile, 'ydot', 'pp', 'u1', 'u2', 'y', 'y1', 'mu', '-mat');
    else
        save(savefile, 'yhat', 'pp', 'u', 'y', 'y1', 'mu', '-mat');
    end
end


Error_plots=0;
Estimate_plots=1;
Par_sigs=0;
%%%%%%%%%%%%%%%%%%%%%%%%%

if Error_plots==1
    figure
    subplot(4,1,1)
%     plot(t,e0(:,i),'-',t*output_amnt,esim(:,i), '--')
    xlabel('Time (s)')
    ylabel('\epsilon')
    subplot(4,1,2)
    plot(t,dele(:,1))
    %plot(t,dele(:,1,3),t,dele(:,1,4),t,dele(:,1,5),t,dele(:,1,6),t,dele(:,1,7),t,dele(:,1,8),t,dele(:,1,9),t,dele(:,1,10))
    xlabel('Time (s)')
    ylabel('E_{L}')
    subplot(4,1,3)
    plot(t,dele(:,2))
    % plot(t,dele(:,2,3),t,dele(:,2,4),t,dele(:,2,5),t,dele(:,2,6),t,dele(:,2,7),t,dele(:,2,8),t,dele(:,2,9),t,dele(:,2,10))
    xlabel('Time (s)')
    ylabel('E_{R_0}')
    subplot(4,1,4)
    plot(t,dele(:,3))
    %     plot(t,dele(:,3,3),t,dele(:,3,4),t,dele(:,3,5),t,dele(:,3,6),t,dele(:,3,7),t,dele(:,3,8),t,dele(:,3,9),t,dele(:,3,10))
    xlabel('Time (s)')
    ylabel('E_{C_2}')
    %     subplot(5,1,5)
    %     plot(t,esim(:,3),t,esim(:,4),t,esim(:,5),t,esim(:,6),t,esim(:,7),t,esim(:,8),t,esim(:,9),t,esim(:,10))
    %     xlabel('Time (sec)')
    %     ylabel('E_k')
end


if Estimate_plots==1
    figure;
   [row,col] = size(pp);
        par_array = zeros(row,col);
               
        
%         for kp=1:length(par_true)
%                  par_array(kp,1:length(pp)) = par_true(kp);
%         end
        
        for kp = 1:row-1
            subplot(row,1,kp);
            plot(pp(kp,:))
            hold on
           % plot(par_array(kp,:),'r')
            ylabel(['Par' num2str(kp)])
           
        end
end

