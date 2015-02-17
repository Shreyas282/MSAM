function mu = GetMu(Y,y_ave,p,DM,beta)
%INPUTS 
%Y: target
%y_ave - candidate and its perturbed states
% start: starting mu
% scalar: makes adaptation of mu more or less aggressive. 2 elements, one
% for the sign accordance and one for exploration/exploitation tradeoff
% DM: the lscov weights calculated for the perturbed cases. 

% updates & discussion
% 7/23/2013: originally the objective of the mu function was to increase
% the confidence size when the perturbed model was more correlated with the
% target than the candidate. However, there are two problems: 1) this did not 
% take into account the fact that the structural error from lscov might be
% negative (the sign of the perturbation), and 2) it encourages large step
% sizes to be taken which may or may not be a good thing. the new approach
% is to have two weighting factors: one that adds confidence (increases mu)
% when the sign of the lscov result matches the change in correlation (i.e.
% increase mu when the step is in the direction of better correlation) and 
% one that uniformly penalizes overly large steps. 

num = length(beta); % number of model compartments
start = p.mod_adapt.mustart;
scalar = p.mod_adapt.agg;
mu_min = p.mod_adapt.mu_min;
if nargin<4
    scalar(1)=100;
    scalar(2)=10;
    if nargin<3
    start=.2;
    end
end
mu = ones(num,1); %separate mu for each term
for x=1:num
    mu(x,:) = mu(x,:)*start(x); % mu starting value
end

adj = zeros(size(y_ave,2)-1,1);


y = y_ave(:,1);
y_pert = y_ave(:,2:end);

%nom_corr = max(xcorr(Y,y,'coeff'));
tmp = corrcoef(Y,y);
nom_corr = tmp(1,2);
nom_err = sum(abs(Y-y));
    
for k=1:size(y_pert,2)
%    pert_corr(k) = max(xcorr(Y,y_pert(:,k),'coeff'));
     tmp = corrcoef(Y,y_pert(:,k));
     pert_corr(k) = tmp(1,2);
     pert_err(k) = sum(abs(Y-y_pert(:,k)));
    
%     adj(k) = (sum(pert_corr(:,k))-sum(nom_corr))/sum(nom_corr);
%     adj(k) = (nom_err-pert_err(k))/nom_err;
    adj(k) = (pert_corr(k)-nom_corr)/nom_corr; 
    ab_adj(k) = abs(adj(k)); % size of perturbances
    

        
    adj2(k) = abs(DM(k)/beta(k)); % ratio of chosen step to perturbance size
%      adj2(k) = log(abs(DM(k)/beta(k))); % ratio of chosen step to perturbance size
%     a(k) = max(pert_corr(:,k));
%     b(k) = max(nom_corr);
%     if adj(k)<0
%         mu(k)=.5;
%         if adj2(k)<0
%             mu(k)=.5*mu(k);
%         end
%     elseif adj2(k)<0
%         mu(k)=.25*mu(k);
%     end
    
%     mu(k) = abs(mu(k)*adj(k)*adj2(k));
%     mu(k) = mu(k)+mu(k)*adj(k)*scalar;
    if sign(adj(k)) == sign(DM(k))
        mu(k) = mu(k)+mu(k)*ab_adj(k)*scalar(1) - mu(k)*adj2(k)*scalar(2);
    else
        mu(k) = mu(k)-mu(k)*ab_adj(k)*scalar(1) - mu(k)*adj2(k)*scalar(2);
    end
        
    if mu(k) < mu_min % do not let mu change the sign of the gamma step; the least squares should handle that
        mu(k) = mu_min;
    elseif mu(k) > 1
        mu(k) = 1;
    end
end

% 
% figure;
% subplot(121)
% plot(Y,'k'); hold on;
% plot(y,'g');
% plot(y_pert(:,1),'--b');
% plot(y_pert(:,2),'--r');
% 
% subplot(122)
% plot(nom_corr,'k'); hold on;
% plot(pert_corr(:,1),'--b');
% plot(pert_corr(:,2),'--r');
% 
% set(gcf,'position',[415 175 1130 420]);
% 
% pause(.001);
