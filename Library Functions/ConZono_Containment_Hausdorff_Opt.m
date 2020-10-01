function [z_s] = ConZono_Containment_Hausdorff_Opt(z,x)

% ConZono_Containment_Hausdorff_Opt computes an inner-approximating zonotope based
% on the pre-defined zonotope structure for a constrained zonotope by minimizing
% hausdorff distance.

% Inputs: x - constrained zonotope in CG-Rep as a struct variable with c_1,G_1,A_1,b_1
% c_1-center of zonotope, G_1-generator matrix,A_1,b_1 satisfying 
% x = {c_1 + G_1\xi_1, ||\xi_1||_{\infty} <= 1, A_1\xi_1 = b_1 (constraints)}
% z - Zonotope structure in CG-Rep as a struct variable with c_2,G_2
% z = {c_2 + G_2\xi_2, ||\xi_2||_{\infty} <= 1}

% Returns a zonotope z_s which is a scaled version of z but with a different center

n_dof = size(x.G,1);

x1.c = z.c;
x1.G = z.G;
x1.ng = size(x1.G,2);
x1.P = [eye(x1.ng);-eye(x1.ng)];
x1.q = ones(2*x1.ng,1);
x1.nq = size(x1.P,1);
x1.nP = size(x1.P,2);

T = null(x.A,'r'); % Computes the right null space of A.
s = pinv(x.A)*x.b;
P = Polyhedron('H',[T ones(size(s,1),1)-s; -T ones(size(s,1),1)+s]);
x2.c = x.c+x.G*s;
x2.G = x.G*T;
x2.ng = size(x2.G,2);
x2.P = P.H(:,1:end-1);
x2.q = P.H(:,end);
x2.nq = size(x2.P,1);
x2.nP = size(x2.P,2);

P_inf = [eye(n_dof);-eye(n_dof)];
q_inf = ones(2*n_dof,1);

phi = sdpvar(x1.ng,1);
cm = sdpvar(n_dof,1);
d = sdpvar(1,1);
Gamma_a = sdpvar(x2.ng,x1.ng,'full');
beta_a = sdpvar(x2.ng,1);
Lambda_a =  sdpvar(x2.nq,x1.nq,'full');

Gamma_b = sdpvar(n_dof,x1.ng,'full');
beta_b = sdpvar(n_dof,1);
Lambda_b =  sdpvar(2*n_dof,x1.nq,'full');

Gamma_c = sdpvar(x1.ng,x2.ng,'full');
beta_c = sdpvar(x1.ng,1);
Lambda_c =  sdpvar(x1.nq,x2.nq,'full');

Gamma_d = sdpvar(n_dof,x2.ng,'full');
beta_d = sdpvar(n_dof,1);
Lambda_d =  sdpvar(2*n_dof,x2.nq,'full');

cons = [];
cons = [cons, 0 <= phi]; % Maximizes the scaling variable
cons = [cons, 0 <= d]; % Minimizes the hausdorff distance
cons = [cons, Lambda_a >= 0];
cons = [cons, Lambda_b >= 0];
cons = [cons, Lambda_c >= 0];
cons = [cons, Lambda_d >= 0];

% Enforces z_s \subset x
cons = [cons, x1.G*diag(phi) == x2.G*Gamma_a + Gamma_b];
cons = [cons, x2.c - cm == x2.G*beta_a + beta_b];
cons = [cons, Lambda_a*x1.P == x2.P*Gamma_a];
cons = [cons, Lambda_a*x1.q <= x2.q + x2.P*beta_a];
cons = [cons, Lambda_b*x1.P == P_inf*Gamma_b];
% cons = [cons, Lambda_b*x1.q <= d*q_inf + P_inf*beta_b];
cons = [cons, Lambda_b*x1.q <= 0*q_inf + P_inf*beta_b];

% % Enforces z_s \oplus d*B \superset x (d - Hausdorff distance, B - Unit hypercube)
cons = [cons, x2.G == x1.G*Gamma_c + Gamma_d];
cons = [cons, cm - x2.c == x1.G*beta_c + beta_d];
cons = [cons, Lambda_c*x2.P == x1.P*Gamma_c];
cons = [cons, Lambda_c*x2.q <= diag([phi;phi])*x1.q + x1.P*beta_c];
cons = [cons, Lambda_d*x2.P == P_inf*Gamma_d];
cons = [cons, Lambda_d*x2.q <= d*q_inf + P_inf*beta_d];

objs = d;
opts = sdpsettings('solver','gurobi','verbose',0);
[sol] = optimize(cons,objs,opts); % Solves the LP

z_s = z;
z_s.c = value(cm);
z_s.G = z.G*diag(value(phi));

end