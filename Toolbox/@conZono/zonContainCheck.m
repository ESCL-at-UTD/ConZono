function cons_out = zonContainCheck(obj1,obj2,cons_in)
% Check if obj1 is contained in obj2

if obj1.n ~= obj2.n
    disp(['Cannot check containment of a AH-polytope in ',num2str(obj1.n),...
        ' dimensions to a AH-polytope in ',num2str(obj2.n),' dimensions!'])
else

	Gamma = sdpvar(obj2.nG,obj1.nG,'full');
	beta = sdpvar(obj2.nG,1);

% 	cons_in = [];
	cons_in = [cons_in, [obj1.G] == obj2.G*Gamma];
	cons_in = [cons_in, obj2.c-(obj1.c) == obj2.G*beta];
	cons_in = [cons_in, abs(Gamma)*ones(obj1.nG,1)+abs(beta) <= ones(obj2.nG,1)];
    cons_out = cons_in;
%    objs = 0; % Feasibility problem
%    opts = sdpsettings('solver','linprog','verbose',0); % Solver settings
%    sol = optimize(cons,objs,opts); % Solves the LP.
%    if sol.problem == 0
%        out = 1;
%    else
%        out = 0;
%    end

end