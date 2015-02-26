function eqn_form = GetEqnForm(models,perts)
%returns equation string from model structure. 
%INPUT
% model: structure
% pert: index of perturbation in model.terms structure array
% old = digits(4);
if nargin<2
    
    eqn_form = vpa(sum([models(1).terms(:).val]),4);
else
    pertoibed = models(1).terms;
   for i=1:length(perts)         
        pertoibed(perts(i)).val = pertoibed(perts(i)).val*models(i+1).terms(perts(i)).pert;
   end
    
    eqn_form = vpa(sum([pertoibed(:).val]),4);
% old notation:
%     eqn_sym = model.t1.val + model.t2.val*model.t2.pert + model.t3.val;
end
% digits(old);
end