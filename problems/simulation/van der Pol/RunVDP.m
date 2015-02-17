function out = RunVDP(model,p)

  
%     u=0;
%     input_set = p.simulation.input_set;% 1:impulse; 2:step; 3:random
    out.x1 = [];
    out.dx1 = [];
    out.ddx1 = [];
    
    time_span=p.simulation.time_span;
    samp_time=time_span/(p.simulation.ndata-1);
    step_size = p.simulation.step_size;
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
    outstr = model{1};
    try
            yout=[];
            for h=1:size(p.simulation.IC,1)
                for i = 1:size(p.simulation.IC,2)
                IC = p.simulation.IC(h,i,:);
                
                sim(p.simulation.sim_model,[],options);
                
                out.x1 = [out.x1;x1];
                
                out.ddx1 = [out.ddx1;ddx1];
                out.dx1 = [out.dx1; dx1];
                end
            end
            y_tmp = yout;
    catch
        disp('Run Model simulation error.');
        keyboard
    end
%     if p.simulation.partitioning==1
%         out.x1 = x1;
%         out.ddx1 = ddx1;
%         out.dx1 = dx1;
% %         out.u_in = reshape(u,[p.simulation.ndata,1]);
% %         out(4) = time;
%     else
%         out = y_tmp;
%     end

%     keyboard
end
