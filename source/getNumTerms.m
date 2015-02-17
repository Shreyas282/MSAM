function number = getNumTerms(eqn_sym)
% splits symbolic equation into linear combination of expanded terms and
% returns the number of terms. 

f = eqn_sym;


fexchar = char(expand(f));

plusminus = regexp(fexchar,{'+','-'});
parenth = regexp(fexchar,{'(',')'});

pm = [plusminus{:}];
pm = sort(pm);
p_in = [parenth{1}];
p_out = parenth{2};

% pms between ) and ( are the spots we can split the equation into terms
for k = 1:length(p_in)
    i1 = pm>p_in(k);
    i2 = pm<p_out(k);
    pm(i1==i2) = 0;
end
pm=pm(pm~=0);
pm=pm(pm~=1);
number = length(pm(pm~=0))+1;
% if length(pm) == 1
%     number = 1;
% end
end