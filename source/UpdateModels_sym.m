function models = UpdateModels_sym(models,g,p)
% update symbolic models with perturbances perts with exponents g

can_mod=models(1);
pert_mods = models(2:end);

opts=find(strcmpi({can_mod.terms(:).type},'int'));
try
    if p.perturb_extvars
        if p.perturb_extvars
            opts = [opts; find(strcmpi({can_mod.terms(:).type},'ext'))];
        end
    end
end
for j = 1:length(pert_mods)
    can_mod.terms(opts(j)).val = pert_mods(j).terms(opts(j)).val*...
                                 pert_mods(j).terms(opts(j)).pert;
    can_mod.terms(opts(j)).val = subs(can_mod.terms(opts(j)).val,...
                                      {['b' num2str(j)]},[g(j)]);
                                  
end

can_mod.eqn_sym = GetEqnSym(can_mod);
can_mod.eqn_str = GetEqnStr_sym(can_mod,p.allvars);


for j = 1:length(pert_mods)
    for count=1:length([can_mod.terms])
        pert_mods(j).terms(count).val = can_mod.terms(count).val; 
    end
    pert_mods(j).eqn_sym = GetEqnSym(pert_mods(j),opts(j));
    pert_mods(j).eqn_str = GetEqnStr_sym(pert_mods(j),p.allvars);
end


models = [can_mod; pert_mods];

end