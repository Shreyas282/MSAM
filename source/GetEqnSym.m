function eqn_sym = GetEqnSym(model,p,pert)
%returns equation string from model structure. 
%INPUT
% model: structure
% pert: index of perturbation in model.terms structure array

if nargin<3
   
    eqn_sym = sum([model.terms(:).val]);
else
            
    pertoibed = model.terms;
    pertoibed(pert).val = pertoibed(pert).val*pertoibed(pert).pert;
    eqn_sym = sum([pertoibed(:).val]);
% old notation:
%     eqn_sym = model.t1.val + model.t2.val*model.t2.pert + model.t3.val;
end
end