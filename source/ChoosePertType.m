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
c = 1;
for count=1:length(pert_index)
    if pert_index(count)~=0
        eval(['syms g' num2str(c)]);

        %choose perturbations randomly
    %      pertnum(count) = intrand(1,length(pert_options));
         %pertnum(count) = pert_index(count);
    % %     debugging
    %         pertnum(1) = 2;
    %         pertnum(2) = 1;
    %         pertnum(3) = 3;

        pertchoice(c) = pert_options(pert_index(count));

        chosen(c) = candidate;

        if p.mod_adapt.useabs   
            chosen(c).terms(count).pert = eval(['abs(pertchoice(c).^g' num2str(c) ');']); 
        else
            chosen(c).terms(count).pert = eval(['pertchoice(c).^g' num2str(c) ';']); 
        end
        chosen(c).terms(count).gamma = 0;

        chosen(c).eqn_sym = GetEqnSym(chosen(c),p);
        chosen(c).eqn_str = GetEqnStr_sym(chosen(c),p.allvars,p);
        chosen(c).eqn_form = chosen(c).eqn_sym;
        c=c+1;
    end
end
c=1;
for count=1:length(pert_index)
    if pert_index(count)~=0
     % update candidate model
        candidate.terms(count).gamma = 0;
        %check if perturbance matches the term
        if ~isempty(strfind(char(candidate.terms(count).val),char(pertchoice(c)))) 
            for k=1:length(p.cons(1,:))
                eval(['syms ' p.cons{1,k}])
            end
            for k=1:length(p.intvars)
                eval(['syms ' char(p.intvars(k))])
            end
    %         tmp = regexprep(char(candidate.terms(term_options(count)).val),...
    %          char(pertchoice(count)),['sign(' char(pertchoice(count)) ...
    %          ')*' char(pertchoice(count)) '^(1+g' num2str(count) ')']);
            tmp = regexprep(char(candidate.terms(count).val),...
                  char(pertchoice(c)),['abs(' char(pertchoice(c)) ')^(1+g' ...
                  num2str(c) ')']);
            candidate.terms(count).val = eval(tmp);
        else
        candidate.terms(count).val = ...
        candidate.terms(count).val*chosen(c).terms(count).pert;
        end
        c=c+1;
    end
end
candidate.eqn_sym = GetEqnSym(candidate,p);
candidate.eqn_str = GetEqnStr_sym(candidate,p.allvars,p);
candidate.eqn_form = GetEqnForm([candidate,chosen]);
chosen=chosen';
% % put gamma terms into candidate
% candidate.eqn_sym = GetEqnSym(candidate);
% candidate.eqn_str = GetEqnStr_sym(candidate);

end