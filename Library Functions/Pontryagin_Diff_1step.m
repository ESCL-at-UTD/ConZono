function [Opt_solve] = Pontryagin_Diff_1step(Zm,Zs)
% Pontryagin_Diff_1step computes an exported optimization model that computes
% the Pontryagin difference between zonotopes.

% Inputs: Z1 - Zonotope in CG-Rep (as a struct variable with c_1,
% G_1) parameters satisfying Z1 = {c_1 +G_1\xi_1,
% ||\xi_1||_{\infty} \leq 1}. 
% Z2 - Zonotope in CG-Rep (as a struct variable with c_2, G_2) parameters 
% satisfying Z2 = {c_2 +G_2\xi_2, ||\xi_2||_{\infty} \leq 1}. 

% Returns a exported optimization model 

n_dof = size(Zm.G,1);

% Zonotope Containment Problem 
clear('yalmip')
cd = sdpvar(n_dof,1);
S = sdpvar(size(Zm.G,2)+size(Zs.G,2),1);
Gamma = sdpvar(size(Zm.G,2),size(Zm.G,2)+2*size(Zs.G,2),'full');
beta = sdpvar(size(Zm.G,2),1);

cons = [];
cons = [cons, [[Zm.G Zs.G]*diag(S) Zs.G] == Zm.G*Gamma];
cons = [cons, Zm.c-(cd+Zs.c) == Zm.G*beta];
cons = [cons, abs(Gamma)*ones(size(Zm.G,2)+2*size(Zs.G,2),1)+abs(beta) <= ones(size(Zm.G,2),1)];
cons = [cons, S >= 0];

obj = norm([Zm.G Zs.G]*(1-S),'inf'); % Norm p = 1, 2, \infty
opts = sdpsettings('solver','gurobi','verbose',0);
Opt_solve = optimizer(cons,obj,opts,[],{cd,S}); % Create exported optimization model object

end