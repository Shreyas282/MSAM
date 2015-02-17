function [chosen,candidate,pertnum] = ChoosePertType(candidate,p,pert_index)

% candidate: model structure
% options: perturbance options (variables considered)
if p.mod_adapt.useabs
    pert_options=p.absintvars;
else
    pert_options=p.intvars;
end

% term options are those that are marked as internal variables
term_options=find(strcmpi({candidate.terms(:).type},'int'));

for count=1:length(term_options)
    eval(['syms b' num2str(count)]);

    %choose perturbations randomly
%      pertnum(count) = intrand(1,length(pert_options));
     pertnum(count) = pert_index(count);
% %     debugging
%         pertnum(1) = 2;
%         pertnum(2) = 1;
%         pertnum(3) = 3;
    
    pertchoice(count) = pert_options(pertnum(count));
    
    chosen(count) = candidate;
       
    chosen(count).terms(term_options(count)).pert = eval(['pertchoice(count).^b' num2str(count) ';']); 
    
    chosen(count).eqn_sym = GetEqnSym(chosen(count),term_options(count));
    chosen(count).eqn_str = GetEqnStr_sym(chosen(count),p.allvars);
    
end

chosen=chosen';
% % put gamma terms into candidate
% candidate.eqn_sym = GetEqnSym(candidate);
% candidate.eqn_str = GetEqnStr_sym(candidate);

end