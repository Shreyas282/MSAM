function eqn_str = GetEqnStr_sym(model,vars,p,opt)
%returns equation string from model structure. 
% if nargin<2 % NOTE: all this functionality has been passed to GetEqnSym.m
%     fchar = char(model.t1.val+model.t2.val+model.t3.val);
%     eqn_str = regexprep(fchar,{'dx1','x1','u_in'},{'u(1)','u(2)','u(3)'});
% elseif pert==1
%     fchar = char(model.t1.val*model.t1.pert+model.t2.val+model.t3.val);
%     eqn_str = regexprep(fchar,{'dx1','x1','u_in'},{'u(1)','u(2)','u(3)'});
% elseif pert==2
%     fchar = char(model.t1.val+model.t2.val*model.t2.pert+model.t3.val);
%     eqn_str = regexprep(fchar,{'dx1','x1','u_in'},{'u(1)','u(2)','u(3)'});
% end
if nargin<4
    opt = 'mod';
end
if strcmpi(opt,'mod')       
    fchar = char(model.eqn_sym);
elseif strcmpi(opt,'eqn')
%     tmp = regexprep(char(candidate.terms(term_options(count)).val),...
%               char(pertchoice(count)),[char(pertchoice(count)) '^(1+g' ...
%               num2str(count) ')']);
    fchar = char(model);
else
    disp(['BAD INPUT TO GETEQNSTR_SYM.M']);
end
% eqn_str = regexprep(fchar,{'dx1','x1','u_in'},{'u(1)','u(2)','u(3)'});
for j=1:length(vars)
    reps(j) = {['u(' num2str(j) ')']};
    charvars(j) = {char(vars(j))};
    charvars_abs(j) = {['sign(' char(vars(j)) ')*abs(' char(vars(j)) ')']};
end

if p.mod_adapt.useabs
    for i = 1:p.num_terms
        for j = 1:length(vars)
            if ~isempty(strfind(fchar,['abs(' char(vars(j)) ')^(g' num2str(i) ' + 1)']))
                fchar = strrep(fchar,['abs(' char(vars(j)) ')^(g' num2str(i) ' + 1)'], ...
                    ['sign(' char(vars(j)) ')*abs(' char(vars(j)) ')^(g' num2str(i) ' + 1)']);
            end
        end
    end
end

eqn_str = regexprep(fchar,charvars,reps);
end