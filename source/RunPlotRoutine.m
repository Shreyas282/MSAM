function perttitle = RunPlotRoutine(h,a,x,its,y_ave,p,colors,...
                            sum_abs_error,max_corr,pertnum,p_t,phi,gamma,...
                            error)
% plot in loop results.
    leggy = {};
    leg={};
    %exp = p.mod_adapt.exp;
    
    
    
    figure(h);

    for count=1:length(y_ave(1,2:end,1))
        leggy{count} = ['$\hat{y}+\delta M_' num2str(count) '$'];
        leg{count} = ['$\gamma_' num2str(count) '$'];
%         leg{2*count} = ['$\gamma_' num2str(count) '$'];
    end                
    leggy = {'$y^*$','$\hat{y}$',leggy{:}};

%             set(gcf,'position',[200,200,1200,600]);

    set(gcf,'units','normalized','position',[.2 .2 .6 .6]);
    
    subplot(2,3,[1 4])
    plot(p.Y,'k'); hold on;
    try
    plot(y_ave(:,1,x),'-.xk','Linewidth',2); 
    catch
        keyboard
    end
    for count=1:length(y_ave(1,2:end,1))
       plot(y_ave(:,count+1,x),'linestyle','--','color',colors(count));

    end
    hold off;
    l=legend(leggy,'Location','Best');
    set(l,'interpreter','latex','Fontsize',14);
    set(gca,'FontSize',14);
    
    subplot(232)
    AX = plotyy(1:x,sum_abs_error,1:x,max_corr);
    pertlist='';
    pertnum(pertnum==0)=[];
    for count=1:length(pertnum(:))
        pertlist=strcat(pertlist,',',char(p.intvars(pertnum(count))));
    end
    perttitle = strcat('Perturbances [', pertlist(2:end), ']');
    title(perttitle,'FontSize',16,'interpreter','latex');
    xlabel('iterations','FontSize',14,'interpreter','latex');
    set(get(AX(1),'Ylabel'),'String','$\sum{|y^*-\hat{y}|}$','interpreter','latex','FontSize',16) 
    set(get(AX(2),'Ylabel'),'String','COV($y^*$,$\hat{y}$)','interpreter','latex','FontSize',16) 
    xlim([0 its+1]);
    ylim([0 max([sum_abs_error;max_corr])]);

    hold on

    % NLS performance
    subplot(235)
    leggy2{1}='$\epsilon(t)$';
    plot(error(:,x),'k'); hold on;
%     if ~p.mod_adapt.MC_NLS_plot
%         for count=1:length(p_t(:))
%             plot(phi(:,count)*p_t(count),'linestyle','--','color',colors(end-count));
%             leggy2{count+1} = ['$\hat{\Delta M} \Phi(t)_' num2str(count) '$'];
%         end
%     end
    leggy2{end+1} = '$\sum{\Delta \hat{M} \Phi(t)}$';
    plot(phi*p_t(:)); 
    if p.mod_adapt.MC_NLS_plot 
%         for count2 = 1:size(p.dydtheta,2)-1
%             for count=1:p.simulation.ndata
%                 dydtheta_norm(count,count2) = norm(...
%                         p.mod_adapt.dydtheta_star(count,:)-...
%                         reshape(p.dydtheta(count,count2+1,:),[1 size(p.dydtheta,3)]));    
%             end
%             dMCstar(count2) = sum(dydtheta_norm(:,count2))/p.mod_adapt.MCstar;
%             
%         end
        dMCstar = (p.true_gamma'-gamma(:,end));
        plot(phi*(p.MC(1)-repmat(p.mod_adapt.MCstar,...
            size(p.MC(1))))/p.mod_adapt.MCstar,'--k','Linewidth',2);
        plot(phi*dMCstar,'--k','Linewidth',2);
        leggy2{end+1} = '$\sum{\Delta M^* \Phi(t)}$';
    end
        hold off;
    l=legend(leggy2{:});%,...
            %'Location','NorthOutside','Orientation','horizontal');
    set(l,'interpreter','latex','Fontsize',14);
    
        

    subplot(2,3,[3 6])        
    for count=1:length(gamma(:,1))
        plot(1:x,gamma(count,1:x),'color',colors(count)); hold on
        % this line plot is not available when gamma and target are different lengths
%         line([0 x+1],[p.mod_adapt.exp(count) p.mod_adapt.exp(count)],'Color',colors(count),'linewidth',1.5,'linestyle','--');
    end
    hold off
    xlabel('iterations','FontSize',14,'interpreter','latex');
    xlim([0 its+1]);
    l=legend(leg);
    set(l,'interpreter','latex','fontsize',16);
%     ylim([0 max([max(gamma)])*2 ]);
    
    pause(.5);
end