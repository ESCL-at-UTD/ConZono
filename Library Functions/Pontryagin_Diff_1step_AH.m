function [Opt_solve] = Pontryagin_Diff_1step_AH(Zm,Zs)

% Pontryagin_Diff_1step_AH computes an exported optimization model that computes
% the Pontryagin difference between zonotopes using AH polytope containment

% Inputs: Zm - Zonotope in CG-Rep (as a struct variable with c_m,
% G_m) parameters satisfying Zm = {c_m +G_m\xi_m,
% ||\xi_m||_{\infty} \leq 1}. 
% Zs - Zonotope in CG-Rep (as a struct variable with c_s, G_s) parameters 
% satisfying Zs = {c_s +G_s\xi_s, ||\xi_s||_{\infty} \leq 1}. 

% Returns a exported optimization model 

n_dof = size(Zm.G,1);
ng_min = size(Zm.G,2);
ng_sub = size(Zs.G,2);

% Zonotope Containment Problem 
clear('yalmip')
cd = sdpvar(n_dof,1);
S = sdpvar(ng_min+ng_sub,1);
Gamma = sdpvar(ng_min,ng_min+ 2*ng_sub,'full');
beta = sdpvar(ng_min,1);
Lambda = sdpvar(2*ng_min,2*ng_min + 4*ng_sub,'full');

H_1 = [eye(ng_min + 2*ng_sub); -eye(ng_min + 2*ng_sub)];
h_1 = [ones(2*ng_min + 4*ng_sub,1)];

H_2 = [eye(ng_min); -eye(ng_min)];
h_2 = [ones(2*ng_min,1)];

cons_fix = [];
cons_fix = [cons_fix, Lambda*H_1 == H_2*Gamma];
cons_fix = [cons_fix, Lambda*h_1 <= h_2 + H_2*beta];
cons_fix = [cons_fix, Lambda >= 0];
cons_fix = [cons_fix, S >= 0];
cons_chg = [];
cons_chg = [cons_chg, [[Zm.G Zs.G]*diag(S) Zs.G] == Zm.G*Gamma];
cons_chg = [cons_chg, Zm.c-(cd+Zs.c) == Zm.G*beta];
cons = [cons_fix,cons_chg];

obj = norm([Zm.G Zs.G]*(1-S),'inf'); % Norm p =1, 2, \infty
opts = sdpsettings('solver','gurobi','verbose',0);
Opt_solve = optimizer(cons,obj,opts,[],{cd,S}); % Create exported optimization model object

end