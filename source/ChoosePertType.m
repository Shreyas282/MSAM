function [chosen,candidate] = ChoosePertType(candidate,p,pert_index)

% candidate: model structure
% options: perturbance options (variables considered)
% if p.mod_adapt.useabs
    % pert_options=p.absintvars;
% else
    % pert_options=p.intvars;
% end
pert_options = p.intvars;
% term options default those that are marked as internal variables
term_options=find(strcmpi({candidate.terms(:).type},'int'));
try
if p.perturb_extvars
    term_options = [term_options; find(strcmpi({candidate.terms(:).type},'ext'))];
end
end
for count=1:length(term_options)
    
    eval(['syms g' num2str(count)]);

    %choose perturbations randomly
%      pertnum(count) = intrand(1,length(pert_options));
     %pertnum(count) = pert_index(count);
% %     debugging
%         pertnum(1) = 2;
%         pertnum(2) = 1;
%         pertnum(3) = 3;
    
    pertchoice(count) = pert_options(pert_index(count));
    
    chosen(count) = candidate;
	
    if p.mod_adapt.useabs   
		chosen(count).terms(term_options(count)).pert = eval(['abs(pertchoice(count).^g' num2str(count) ');']); 
    else
		chosen(count).terms(term_options(count)).pert = eval(['pertchoice(count).^g' num2str(count) ';']); 
    end
	chosen(count).terms(term_options(count)).gamma = 0;
    
    chosen(count).eqn_sym = GetEqnSym(chosen(count),term_options(count));
    chosen(count).eqn_str = GetEqnStr_sym(chosen(count),p.allvars,p);
    chosen(count).eqn_form = chosen(count).eqn_sym;
end
for count=1:length(term_options)
 % update candidate model
    candidate.terms(term_options(count)).gamma = 0;
    %check if perturbance matches the term
    if ~isempty(strfind(char(candidate.terms(term_options(count)).val),char(pertchoice(count)))) 
        for k=1:length(p.cons(1,:))
            eval(['syms ' p.cons{1,k}])
        end
        for k=1:length(p.intvars)
            eval(['syms ' char(p.intvars(k))])
        end
%         tmp = regexprep(char(candidate.terms(term_options(count)).val),...
%          char(pertchoice(count)),['sign(' char(pertchoice(count)) ...
%          ')*' char(pertchoice(count)) '^(1+g' num2str(count) ')']);
        tmp = regexprep(char(candidate.terms(term_options(count)).val),...
              char(pertchoice(count)),['abs(' char(pertchoice(count)) ')^(1+g' ...
              num2str(count) ')']);
        candidate.terms(term_options(count)).val = eval(tmp);
    else
    candidate.terms(term_options(count)).val = ...
    candidate.terms(term_options(count)).val*chosen(count).terms(term_options(count)).pert;
    end
end
candidate.eqn_sym = GetEqnSym(candidate,p);
candidate.eqn_str = GetEqnStr_sym(candidate,p.allvars,p);
candidate.eqn_form = GetEqnForm([candidate,chosen],term_options);
chosen=chosen';
% % put gamma terms into candidate
% candidate.eqn_sym = GetEqnSym(candidate);
% candidate.eqn_str = GetEqnStr_sym(candidate);

end