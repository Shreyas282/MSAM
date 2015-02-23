function out = EvalModel(model,p)
% evaluate a model algebraically.
    % OUTPUT
    % out : time series output
    
%   if nargin<3
%       target=0;
%   end
    
  % declare parameter values from stored values
     for count=1:length([p.cons(1,:)])
        eval([p.cons{1,count} ' = ' num2str(p.cons{2,count}) ';']);
     end
     % declare variable values
    for count = 1:length(p.allvars)
        eval([char(p.allvars(count)) ' = p.' char(p.allvars(count)) ';']);
    end
    
    model_str = regexprep(char(model),'\^','\.\^');
    model_str = regexprep(model_str,'\*','\.\*');
    model_str = regexprep(model_str,'\/','\.\/');

    try

       out = eval([model_str ';']);
        pass=1;
    catch

        fprintf(['Model form ' model_str ...
            ' failed to output a valid response in simulink.\n']);
        keyboard
        pass=0;
    end


%     keyboard
end
