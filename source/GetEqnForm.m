function eqn_form = GetEqnForm(models,perts)
%returns equation string from model structure. 
%INPUT
% model: structure
% pert: index of perturbation in model.terms structure array

if nargin<2
    
    eqn_form = sum([models(1).terms(:).val]);
else
    pertoibed = models(1).terms;
   for i=1:length(perts)         
        pertoibed(perts(i)).val = pertoibed(perts(i)).val*models(i+1).terms(perts(i)).pert;
    end
    eqn_form = sum([pertoibed(:).val]);
% old notation:
%     eqn_sym = model.t1.val + model.t2.val*model.t2.pert + model.t3.val;
end
end