function out = RunModel(model,p)

  
    u_in=p.simulation.u_in;
    input_set = p.simulation.input_set;% 1:impulse; 2:step; 3:random
    time_span=p.simulation.time_span;
    samp_time=p.simulation.samp_time;
    IC = p.simulation.IC;
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
    outstr = model{1};
    try
            yout=[];
            if strcmpi(p.simulation.sim_model,'forced_msd_uin')
                for h=1:length(p.simulation.step_size)
                    step_size = p.simulation.step_size(h);
                    sim(p.simulation.sim_model,[],options);
                    yout = [yout;y_tmp];
                end
            else
                sim(p.simulation.sim_model,[],options);
            end
            yout = y_tmp;
    catch
        disp('Run Model simulation error.');
        keyboard
    end
%     if p.simulation.partitioning==1
    out.target = y_tmp;
    if exist('t','var'), out.t = t_tmp; end
    if exist('ddx1','var'),out.ddx1 = ddx1; end
    if exist('dx1','var'), out.dx1 = dx1; end
    if exist('x1','var'), out.x1 = y_tmp; end
    if exist('u_in','var'), out.u_in = u_in; end
%         out(4) = time;
%     else
%         out = y_tmp;
%     end

%     keyboard
end
