function [yhat, effort] = RunCtrlModel_sym(model,p)

% outstr = sim_model;
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
  % define gamma terms
  for i = 1:length(model.terms)
      eval(['g' num2str(i) ' = model.terms(i).gamma;']);
  end
  
outstr = model.eqn_str;  

sim(p.simulation.sim_model,[],options);

yhat = y_tmp;
effort = u;
end