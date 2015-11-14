function [yhat, effort] = RunCtrlModel(sim_model,p)

outstr = sim_model;
time_span = p.simulation.time_span;
samp_time = p.simulation.samp_time;
 options=simset('SrcWorkspace','current','DstWorkspace','current');
 noise_on=0;
 
  if ~isempty(p.cons)
        for count=1:length([p.cons(1,:)])
            if ~isempty(p.cons{1,count})
                eval([p.cons{1,count} ' = ' num2str(p.cons{2,count}) ';']);
            end
        end
  end
    
sim(p.simulation.sim_model,[],options);

yhat = y_tmp;
effort = u;
end