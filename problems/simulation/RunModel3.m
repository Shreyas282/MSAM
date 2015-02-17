function p = RunModel3(models,p)

  
    u=0;
    ndata=p.simulation.ndata;
%     input_set = p.simulation.input_set;% 1:impulse; 2:step; 3:random
    time_span=p.simulation.time_span;
    samp_time=time_span/(p.simulation.ndata-1);
    IC = p.simulation.IC;
    ICs = p.simulation.ICs;
    p.dx1 = [];
    p.x1 =  [];
    p.dx2 = [];
    p.x2 =  [];
    p.dx3 = [];
    p.x3 =  [];
   if strcmpi(p.simulation.sim_model,'forced_msd_uin')
       simin(:,1) = 0:samp_time:time_span; 
       simin(:,2)= cos(1*simin(:,1));
%        [0:249,zeros(1,500),100*ones(1,250)]'+wgn(p.simulation.ndata,1,1);
   end
%     par = [p.cons(2,:)];
    % declare constant values
    if ~isempty(p.cons)
        for count=1:length([p.cons(1,:)])
            if ~isempty(p.cons{1,count})
                eval([p.cons{1,count} ' = ' num2str(p.cons{2,count}) ';']);
            end
        end
    end
    options=simset('SrcWorkspace','current','DstWorkspace','current');
    % set target
%     m=375; mstar =m;
%     c=9800; cstar=c;
%     k=130000; kstar=k;

        outstr1 = models{1};
        outstr2 = models{2};
        outstr3 = models{3};
            yout=[];
            pass=0;
            tries=0;
%    while ~pass    && tries<1000
%        IC = [10*(rand()) 10*(rand()) 10*(rand())];
%        disp(['IC: ' num2str(IC(1)) ' ' num2str(IC(2)) ' ' num2str(IC(3))]);
%      try       
         for i=1:size(ICs,1)
             IC = ICs(i,:);
            sim(p.simulation.sim_model,[],options);
            t=time;
            [FO,goodness] = fit(t,x1,'cubicinterp');
            fx = differentiate(FO,t);
            p.dx1 = [p.dx1; fx];
            p.x1 = [p.x1; x1];
            
            [FO,goodness] = fit(t,x2,'cubicinterp');
            fx = differentiate(FO,t);
            p.dx2 = [p.dx2; fx];
            p.x2 = [p.x2; x2];
            
            [FO,goodness] = fit(t,x3,'cubicinterp');
            fx = differentiate(FO,t);
            p.dx3 = [p.dx3; fx];
            p.x3 = [p.x3; x3];
            
%             plot(t,dx1,'--r'); hold on;
%             plot(t,p.dx1((ndata*(i-1)+1):end),'-.r');
%             plot(t,dx2,'--g');
%             plot(t,p.dx2((ndata*(i-1)+1):end),'-.g');
%             plot(t,dx3,'--k');
%             plot(t,p.dx3((ndata*(i-1)+1):end),'-.k');
% pause(.01)
%             out.dx1 = [out.dx1; dx1];
%             out.x1 =  [out.x1; x1];
%             out.dx2 = [out.dx2; dx2];
%             out.x2 =  [out.x2; x2];
%             out.dx3 = [out.dx3; dx3];
%             out.x3 =  [out.x3; x3];
         end
            pass=1;
%      catch
%         tries=tries+1;
% %         disp('Run Model simulation error.');
% %         keyboard
%      end
%    end
    if p.simulation.partitioning==1
        switch p.simulation.targetstate
            case 1,
                p.Y = p.dx1;
            case 2,
                p.Y = p.dx2;
            case 3,
                p.Y = p.dx3;    
            otherwise,
                disp('Run Model target state error.');
        end
        
    else
        p = y_tmp;
    end

%     keyboard
end
