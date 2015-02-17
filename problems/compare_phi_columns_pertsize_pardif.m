clear
load('Data_snapshot/MSD_NLS_randcount1_mu100_dev50_its10-2_pardif_used_pertsize.mat');

% load('Data_snapshot/MSD_NLS_randcount1_mu100_dev50_its10-2_pardif.mat');



p_true = [1/375; 9800/375;  130000/375];
dT= (p_true-[1/p.cons{2,1}; p.cons{2,2}/p.cons{2,1};p.cons{2,3}/p.cons{2,1}]);
% dT= (-p_true+[1/p.cons{2,1}; p.cons{2,2}/p.cons{2,1};p.cons{2,3}/p.cons{2,1}]);
dT = dT./p_true;
dydt = reshape(dydtheta(:,1,:,:,1),[128 3 its-1]);
% if plot1
%     figure(1); 
%     l=legend('Y','$\hat{y}','$\epsilon_{\Theta}$','$\epsilon_{M}$','$Y-\hat{y} - \epsilon_{\Theta}$');
%     set(l,'interpreter','latex');
% end

true_G = [2;3];
% true_G =[3;2];
for i = 1:size(dydt,3)
    
    phi_CN(i) = 1./phi_index(i+1,1);
   
    c1_norm(i) = norm(phi(:,1,i+1));
     c2_norm(i) = norm(phi(:,2,i+1));
   
 
end
 figure(1); hold on;
 plot(c1_norm,'.r');
 plot(c2_norm,'.b');