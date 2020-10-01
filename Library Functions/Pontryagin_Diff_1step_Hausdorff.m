function [Opt_solve] = Pontryagin_Diff_1step_Hausdorff(Zm,Zs)

% Pontryagin_Diff_1step_Hausdorff computes an exported optimization model that computes
% the Pontryagin difference between zonotopes using zonotope containment by minimizing 
% Hausdorff distance

% Inputs: Zm - Zonotope in CG-Rep (as a struct variable with c_m,
% G_m) parameters satisfying Zm = {c_m +G_m\xi_m,
% ||\xi_m||_{\infty} \leq 1}. 
% Zs - Zonotope in CG-Rep (as a struct variable with c_s, G_s) parameters 
% satisfying Zs = {c_s +G_s\xi_s, ||\xi_s||_{\infty} \leq 1}. 

% Returns a exported optimization model 

n_dof = size(Zm.G,1);
% Zonotope Containment Problem 
clear('yalmip')
cd = sdpvar(n_dof,1);
d = sdpvar(1,1);
S = sdpvar(size(Zm.G,2)+size(Zs.G,2),1);
Gamma = sdpvar(size(Zm.G,2),size(Zm.G,2)+2*size(Zs.G,2),'full');
beta = sdpvar(size(Zm.G,2),1);
Gamma2a = sdpvar(size(Zm.G,2)+size(Zs.G,2),size(Zm.G,2),'full');
Gamma2b = sdpvar(n_dof,size(Zm.G,2),'full');
beta2a = sdpvar(size(Zm.G,2)+size(Zs.G,2),1);
beta2b = sdpvar(n_dof,1);

cons = [];
cons = [cons, [[Zm.G Zs.G]*diag(S) Zs.G] == Zm.G*Gamma];
cons = [cons, Zm.c-(cd+Zs.c) == Zm.G*beta];
cons = [cons, abs(Gamma)*ones(size(Zm.G,2)+2*size(Zs.G,2),1)+abs(beta) <= ones(size(Zm.G,2),1)];

cons = [cons, Zm.G == [Zm.G Zs.G]*Gamma2a + Gamma2b];
cons = [cons, cd - Zm.c == [Zm.G Zs.G]*beta2a + beta2b];
cons = [cons, abs(Gamma2a)*ones(size(Zm.G,2),1) + abs(beta2a) <= S];
cons = [cons, abs(Gamma2b)*ones(size(Zm.G,2),1) + abs(beta2b) <= d];

obj = d;
opts = sdpsettings('solver','gurobi','verbose',0);
Opt_solve = optimizer(cons,obj,opts,[],{cd,S}); % Create exported optimization model object

end