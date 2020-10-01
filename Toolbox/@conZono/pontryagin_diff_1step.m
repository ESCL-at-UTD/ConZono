function out = pontryagin_diff_1step(obj1,obj2)

out = copy(obj1);
obj1.getDimensions
obj2.getDimensions

if obj1.n ~= obj2.n
	disp(['Cannot compute the Pontryagin difference between a zonotope in ',num2str(obj1.n),...
        ' dimensions and a zonotope in ',num2str(obj2.n),' dimensions!'])
else
    cd = sdpvar(obj1.n,1);
    S = sdpvar(obj1.nG+obj2.nG,1);
    Gamma = sdpvar(obj1.nG,obj1.nG+2*obj2.nG,'full');
    beta = sdpvar(obj1.nG,1);
    cons = [];
    cons = [cons, [[obj1.G obj2.G]*diag(S) obj2.G] == obj1.G*Gamma];
    cons = [cons, obj1.c - (cd+obj2.c) == obj1.G*beta];
    cons = [cons, abs(Gamma)*ones(obj1.nG+2*obj2.nG,1) + abs(beta) <= ones(obj1.nG,1)];
    cons = [cons, S >= 0];
    
    obj = norm([obj1.G obj2.G]*(1-S),'inf'); % Norm p = 1, 2, \infty
    opts = sdpsettings('solver','gurobi','verbose',0);
    Opt_solve = optimizer(cons,obj,opts,[],{cd,S}); % Create exported optimization model object
    [OUT_Opt,diagnostics] = Opt_solve([]); % Solves the LP
    out.c = cell2mat(OUT_Opt(1));
    S_val = cell2mat(OUT_Opt(2));
    out.G = [obj1.G obj2.G]*diag(S_val);
%     cd = sdpvar(obj1.n,1);
% 	phi_ = sdpvar(obj1.nG+obj2.nG,1);
% 	obj0 = copy(obj1);
% 	obj0.c = cd + obj2.c;
% 	obj0.G = [obj1.G obj2.G]*diag(phi_);
% 
% 	cons = [phi_ >= 0];
%     
% 	opts = sdpsettings('solver','gurobi','verbose',0);
% 	objs = norm([obj1.G obj2.G]*(1-phi_),'inf');
% 	[sol] = optimize(cons,objs,opts);
% 	out.c = value(cd) + obj2.c;
% 	out.G = [obj1.G obj2.G]*diag(value(phi_));
end
end