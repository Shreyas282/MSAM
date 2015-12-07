function output = getTerms(mod,type,p)
% splits a model's symbolic equation into linear combination of expanded 
% terms and returns a nested structure with terms defined.
% INPUTS
% mod: either a model struct or a symbolic equation
% type: specifies whether input and output are a model struc or equation
%       'mod': model struct input and output. output is struct with term
%       variables
%       'eqn': symbolic equation. output is array of terms. 
% OUTPUT
% type 'mod':
% output.terms : structure array
% output.terms.val : symbolic equation value
% output.terms.pert: empty
% output.terms.type: 'int' if it contains internal variables (like x, dx),
%                    'ext' if it contains external variables (input like u)

chextvars={''};

if strcmpi(type,'mod')
    output=mod;
    f = mod.eqn_sym;
    for count=1:length(p.intvars)
        chintvars(count) = {char(p.intvars(count))};
    end
    for count=1:length(p.extvars)
        chextvars(count) = {char(p.extvars(count))};
    end
elseif strcmpi(type,'eqn')
    f = mod;
else
    fprintf(['Error in getTerms.m']);
    keyboard
end
% declare variables in equation
vars = symvar(f);
for count=1:length(vars)
    eval(['syms ' char(vars(count))]);
end
% syms a b c x y z
fex = expand(f);
fexchar = char(fex);

plusminus = regexp(fexchar,{'+','-'});
parenth = regexp(fexchar,{'(',')'});

pm = [plusminus{:}];
pm = sort(pm);
p_in = [parenth{1}];
p_out = parenth{2};

% pms between ) and ( are the spots we can split the equation into terms
for count = 1:length(p_in)
    i1 = pm>p_in(count);
    i2 = pm<p_out(count);
    pm(i1==i2) = 0;
end

splits = [pm(pm~=0) length(fexchar)+1];
if splits(1)~=1
    splits = [1 splits];
end

if strcmpi(type,'mod')
    
    for j = 1:length(splits)-1
        tmp=0;
        output.terms(j).val = eval(fexchar(splits(j):splits(j+1)-1));
        tmp=regexp([char(output.terms(j).val)],chextvars,'once');
        if sum(tmp{:})==0
            tmp=regexp([char(output.terms(j).val)],chintvars,'once');
            if sum([tmp{:}])==0
                output.terms(j).type = 'cons';
            else
                output.terms(j).type ='int';
            end
        else
            output.terms(j).type = 'ext';
        end
%         output.terms(j) = eval(fexchar(splits(j):splits(j+1)-1));
        output.terms(j).gamma = 0;
    end
elseif strcmpi(type,'eqn')
    for j = 1:length(splits)-1
        output(j) = eval(fexchar(splits(j):splits(j+1)-1));
    end
end
% keyboard
end

