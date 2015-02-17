function [DM,Phi] = StructError(y,y_pert,Y,dMC,opt)
%INPUTS
% y: yhat
% ypert: yhat(M+dM)
% dMC: change in model complexity from y to ypert
%OUTPUT
% DM: lscov of phi and error

if nargin<5
    opt = 'error';
    type=1;
else
    if strcmpi(opt,'error')
        type=1;
    elseif strcmpi(opt,'Ycorr')
        type=2;
    elseif strcmpi(opt,'errcorr')
        type=3;
    elseif strcmpi(opt,'combo')
        type=4;
    else
        fprintf('Input Error (StructError)');
        keyboard
    end
end

if type ~= 2
    Error = Y - y;
end

for k=1:size(y_pert,2)
    Phi(:,k) = (y_pert(:,k) - y)/dMC(k);
%     Phi(:,k) = (y_pert(:,k) - y)/sum(abs(y))/dMC(k);
    % Phi(:,k) = dely(:,k);
    if type == 2
        Y_Phi_corr(:,k) = xcorr(Y,Phi(:,k));
    elseif type==3 || type == 4
        Error_Phi_corr(:,k) = xcorr(Error,Phi(:,k));
    end
end

% Phi_t(:,:,r) = [dely(:,1,r),dely(:,2,r)];
% Phi(:,:) = (ypert(:,:) - repmat(y,1,size(ypert,2)))./dMC(:);
if type == 1 || type == 4

    DM = lscov(Phi,Error);
%     DM = lscov(Phi,Error./dMC);l
    if type == 4
        DM1 = DM;
    end
elseif type == 2
    star_corr = xcorr(Y,Y);
    Y_y_corr =  xcorr(Y,y);
    Error_Y_corr = star_corr-Y_y_corr;
        
    DM = lscov(Y_Phi_corr,Error_Y_corr);
end
if type == 3 || type == 4
    Error_corr = xcorr(Error,Error);
    
    DM = lscov(Error_Phi_corr,Error_corr);
    if type == 4
        DM2 = DM;
    end
end
if type == 4
    DM = mean([DM1, DM2],2);
end
    
% figure(25);
% plot(Error,'k'); hold on;
% for count=1:length(DM)
%     plot(Phi(:,count)*DM(count))
% end
% plot(Phi*DM); hold off;
pause(.1);

end