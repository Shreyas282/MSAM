% This is an iterative least squares adaptation routine

 function cons = par_estimate_VIV(output_focus,yhat)
%  close all
%  clear
% Select the adaptation method; 1=time-based, 2=wavelet-based,
% 3=signature-based

% clear
if nargin<2
    yhat = {'1/A*(QDDOT + eps*(2*pi*St*U/D)*(Q^2-1)*QDOT+(2*pi*St*U/D)^2*Q)'};
end
% yhat = {'1/A*(QDDOT/abs(Q)^(0.0909) + (4 * pi^2 * Q * St^2*U^2*abs(QDDOT)^(6.6580))/(D^2) - (2*pi*QDOT*St*U*eps)/(D*abs(Q)^(0.0277)) + (2*pi*Q^2*QDOT*St*U*eps)/(D*abs(QDDOT)^(0.0026)))'};
% yhat = {'-mu1*(x1)^(2)*dx1+mu2*dx1-mu3*x1'};

model_str = regexprep(yhat{1},'\^','\.\^');
model_str = regexprep(model_str,'\*','\.\*');
model_str = regexprep(model_str,'\/','\.\/');

reload=0;
enableBreak = 1;
Continuous_plot=0;
select_sig=1;
if select_sig==1
   sig_amnt=20; 
end

limit1=1;
limit2=72;
normalize_sensitivity=0; % 1 turns it on, and 0 leaves it off
% mu=0.50;
mu=1;
max_iteration=20;


cascade=1;


Isolate_plots=0;
make_movie=0;
iPlot=0;
cone_remove=1;
%noise_on=1;  %not applied here (must modify model
noise_amp=0.0000005;
Noise_consider=0;
noise_rem=0;

if nargin<1, output_focus=[1]; end
output_amnt=length(output_focus);
outputs=output_amnt;
out_scl=1;
output_select=0;
% isolate_mehod=1;
weight_method='scales';
windowSwitch =0; % enable Hanning windowing of data to squelch spectral leakage
if (windowSwitch == 1)
    hanWindow = hanning(n);
end

out = load('VIV_data.mat');
p.cons = {'A', 'eps'; 
          12,  0.3};

St = out.St;
D = out.D;
% A= out.A;
A = 12;
% eps = out.eps;
eps = 0.3;
y0 = -out.ddy_50s(output_focus,:)'/D;
QDDOT = out.ddforce_50s';
QDOT = out.dforce_50s';
Q = out.force_50s';
U = repmat(out.U',size(Q,1),1);
%% parameters to adapt
par_num = [1:2];
par_start=[p.cons{2,1},p.cons{2,2}];
%  par_start = [mu11_start,mu21_start,mu31_start,gamma(1),gamma(2),gamma(3)];
% par_start = [m1_start,mu21_start,mu31_start,1.7259,2.4250];
 par = par_start;
% for ip = par_num
%     par(ip) = par_start(ip);
% end    

j=length(par_num);
i=1;
% loadfile = ['Forced_nonlinmsd_record.mat'];
% if reload==1
%     load(loadfile);
% end

% esim=zeros(p.simulation.ndata,length(output_focus),max_iteration);
for ii=1:max_iteration
    
    %parameters(:,incremental_adapt)=par_in';
    pp(1:length(par_start),i) = par';   
%     DTheta(:,i)=(par_true(1:length(par))'-pp(1:length(par),i))';
    
    p.cons = {'A', 'eps'; 
           par(1),  par(2)};
    for z=1:size(p.cons,2)
        eval([p.cons{1,z} '= p.cons{2,z};']);
    end
%     A = par(1);
%     eps = par(2);

    y_tmp = (eval([model_str ';']));
    inc=1;
    for kk = output_focus;
        y1(:,inc) = y_tmp(:,kk);
        corr(i,inc) = R2(y0(:,inc),y1(:,inc));
        inc=inc+1;
    end
    
    e0(:,:,i)=y0-y1; % this is the simulation error
    mean_error(i,:) = mean(abs(y0-y1));
    
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
        A = dpar(1,l);
        eps = dpar(2,l);
        
        y_tmp = (eval([model_str ';']));
        
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
    A=p.cons{2,1}; eps=p.cons{2,2}; 
    
    dele1=dele;
    
    
    
    
%     esim(:,:,i)=zeros(size(y0,1),length(output_focus));
    
%     for k = 1:length(output_focus);
%         for l=1:j
% %             del_par(l) = (par_true(par_num(l))-par(par_num(l)))/par_true(par_num(l));
%        %     del_par(l) = (par_true(par_num(l))-par(par_num(l)));
%             esim(:,k,i) = esim(:,k,i) + del_par(l)*dele1(:,k,l);
%         end
%     end
    
      e0Plot = 0;
    
    if (e0Plot == 1)
        figure(2);
               
%        subplot(121); hold off;
%       if ii==1
hold off;
          plot(reshape(y0,size(y0,1)*size(y0,2),1),'b');  
%       end
      hold on;
       plot(reshape(y1,size(y1,1)*size(y1,2),1),'--r');
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
    
    if cascade==1
        
        e0_vect = mat_2_vect(e0(:,:,i));
        e0_vect = e0_vect';
        for l=1:j
            dele2(:,l)=mat_2_vect(dele1(:,:,l));
        end
        outputs=1;
    end
    
    %Gauss-Newton
    if cascade==1
        
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
    end
    

    for l=1:j
        par(par_num(l))=par(par_num(l))+mu*(par_t(l));
        par_trecord(l,i) = par_t(l);
    end
    
    sigDisplay = 3;
    if make_movie ==1
        if i == 1
            blah1 = max(max(coef_mat(:, :,sigDisplay)));
            blah2 = min(min(coef_mat(:, :,sigDisplay)));
            
        end
        figure(1)
        mesh((fliplr(sig(:,:,2).*coef_mat(:, :, sigDisplay)))')
        axis([1 n 1 limit2 blah2 blah1])
        parm_frame(i) = getframe
    end
    
    tmpPar = par';
    i = i+1;

%     save(loadfile, 'i', 'p', 'par', 'par_true', 'pp', 'e0', 'ct_record',  'par_trecord', 'mu', '-mat');

    
    par_error = 0.;
%     for ip=1:j
%         par_error =par_error+ ((par_true(par_num(ip)) - par(ip))/par(par_num(ip)))^2;
%     end
   
    disp(['mean error: ' num2str(mean(mean_error(end,:))) ', R^2 = ' num2str(mean(corr(end,:)))]);
    
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


saveMe =0;

if saveMe ==1

    savefile = ['Results\Forced_nonlinmsd_GN_wBreak_', 'outSelect_',  num2str(output_focus), '_mu_' num2str(mu), 'iter', num2str(max_iteration), '.mat'];
    
end

if saveMe ==1
    

        save(savefile, 'i', 'p', 'par', 'par_true', 'cascade', 'pp', 'e0', 'ct_record',  'par_trecord', 'mu', '-mat');

end


Error_plots=0;
Estimate_plots=0;
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
        
        for kp = 1:row
            subplot(row,1,kp);
            plot(pp(kp,:))
            hold on
            plot(par_array(kp,:),'r')
            ylabel(['Par' num2str(kp)])
           
        end
end

cons = p.cons;
