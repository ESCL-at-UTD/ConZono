function [x_s] = ConZono_Containment_Opt(x_r,x)

% ConZono_Containment_Opt computes an inner-approximating constrained zonotope based
% on the pre-defined constrained zonotope structure for a constrained zonotope by 
% maximizing scaling variables.

% Inputs: x - constrained zonotope in CG-Rep as a struct variable with c_1,G_1,A_1,b_1
% c_1-center of zonotope, G_1-generator matrix, A_1, b_1 satisfying 
% x = {c_1 +G_1\xi_1, ||\xi_1||_{\infty} <= 1, A_1\xi_1 = b_1 (constraints)}
% x_r - Constrained Zonotope structure in CG-Rep with 1 generator and 1 constraint less
% as a struct variable with c_2,G_2, A_2, b_2 satisfying 
% x_r = {c_2 +G_2\xi_2, ||\xi_2||_{\infty} <= 1, A_2\xi_2 = b_2 (constraints)}

% Returns a constrained zonotope x_s which is a scaled version of x_r 
% but with a different center

n_dof = size(x.G,1);

% Affine H-Rep
T = null(x.A,'r'); % Computes the right null space of A.
s = pinv(x.A)*x.b; 
P = Polyhedron('H',[T ones(size(s,1),1)-s; -T ones(size(s,1),1)+s]);

X_Ah.c = x.c+x.G*s;
X_Ah.G = x.G*T;
X_Ah.H = P.H(:,1:end-1);
X_Ah.h = P.H(:,end);
q_x = size(X_Ah.H,1);
n_x = size(X_Ah.H,2);

% Affine H-Rep
T2 = null(x_r.A,'r'); % Computes the right null space of A.
s2 = pinv(x_r.A)*x_r.b;
if isempty(s2)
    s2 = zeros(size(x_r.G,2),1);
end
P2 = Polyhedron('H',[T2 ones(size(s2,1),1)-s2; -T2 ones(size(s2,1),1)+s2]);

X_red_Ah.c = x_r.c+x_r.G*s2;
X_red_Ah.G = x_r.G*T2;
X_red_Ah.H = P2.H(:,1:end-1);
X_red_Ah.h = P2.H(:,end);
q_x_red = size(X_red_Ah.H,1);
n_x_red = size(X_red_Ah.H,2);


yalmip('clear')
Gamma = sdpvar(n_x,n_x_red,'full');
beta = sdpvar(n_x,1);
Lambda =  sdpvar(q_x,q_x_red,'full');
phi = sdpvar(size(x_r.G,2),1);
cm = sdpvar(n_dof,1);

cons = [];
cons = [cons, Lambda >= 0];
cons = [cons, x_r.G*diag(phi)*T2 == X_Ah.G*Gamma];
cons = [cons, X_Ah.c - (cm+x_r.G*diag(phi)*s2) == X_Ah.G*beta];
cons = [cons, Lambda*X_red_Ah.H == X_Ah.H*Gamma];
cons = [cons, Lambda*X_red_Ah.h <= X_Ah.h + X_Ah.H*beta];
cons = [cons, 0 <= phi];

objs = 0;
% objs = norm(1e2-phi,1);
objs = objs + norm(1e2-phi,'inf');
opts = sdpsettings('solver','gurobi','verbose',0);
% opts = sdpsettings('solver','linprog');
[sol] = optimize(cons,objs,opts);

x_s = x_r;
x_s.c = value(cm);
x_s.G = x_r.G*diag(value(phi));

end

