function eqn_str = GetEqnStr_sym(model,vars,opt)
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
if nargin<3
    opt = 'mod';
end
if strcmpi(opt,'mod')
    fchar = char(model.eqn_sym);
elseif strcmpi(opt,'eqn')
    fchar = char(model);
else
    disp(['BAD INPUT TO GETEQNSTR_SYM.M']);
end
% eqn_str = regexprep(fchar,{'dx1','x1','u_in'},{'u(1)','u(2)','u(3)'});
for j=1:length(vars)
    reps(j) = {['u(' num2str(j) ')']};
    charvars(j) = {char(vars(j))};
end
eqn_str = regexprep(fchar,charvars,reps);
% eqn_str = regexprep(fchar,{'x1','x2'},{'u(1)','u(2)'});

% need to generalize this function!!
end