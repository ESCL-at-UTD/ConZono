function out = innerApproxPointContain(obj1,obj2,method,point)
% Method options: 'Hausdorff' or 'Infinity Norm'

if obj1.n ~= obj2.n
    disp(['Cannot approximate a zonotope in ',num2str(obj2.n),...
        ' dimensions by a zonotope in ',num2str(obj1.n),' dimensions!'])
else
    if obj1.nC == 0
        type = 'zono';
    else
        type = 'conZono';
    end
    if isempty(method) && strcmp(type,'zono')
        method = 'Hausdorff';
    elseif isempty(method) && strcmp(type,'conZono')
        method = 'Infinity Norm';
    else
%         method = varargin;
        if  strcmp(method,'Hausdorff') && strcmp(type,'conZono')
            disp('Can only use Hausdorff when inner-approximation is an unconstrained zonotope! Using Infinity Norm instead.')
            method = 'Infinity Norm';
        end
    end
    out = copy(obj1);
    phi_ = sdpvar(obj1.nG,1);
    center_ = sdpvar(obj1.n,1);
    
    obj0 = copy(obj1);
    obj0.c = center_;
    obj0.G = obj1.G*diag(phi_);
    h0 = conZono2AHPoly(obj0);
    h1 = conZono2AHPoly(obj1);
    h2 = conZono2AHPoly(obj2);
    
    cons = [phi_ >= 0];
    opts = sdpsettings('solver','gurobi','verbose',0);
    
    if strcmp(method,'Infinity Norm')
        cons = conContainCheck(h0,h2,cons);
        cons = conPointContain(point,h0,cons); % To add function
        objs = norm(1e2-phi_,'inf');  %%%%%% Hard coded upper-bound
        
    elseif strcmp(method,'Hausdorff')
        d_ = sdpvar(1,1);
        B = conZono;
        B.c = zeros(obj1.n,1);
        B.G = eye(obj1.n);
        B_ah = conZono2AHPoly(B);
        B_d = B_ah*d_;
        
%         h1.c = h0.c;
        h1.c = center_;
        h1.b = diag([phi_;phi_])*h1.b;
        
        cons = [cons, d_ >= 0];
        cons = conContainCheck(h0,h2,cons);
		cons = conPointContain(point,h0,cons);
        cons = conContainCheck(h2,h1+B_d,cons);
        objs = d_;  %% Distance alone did not seem to maximize sets.
        objs = objs + norm(1e2-phi_,'inf');  %%%%%% Hard coded upper-bound
    end
    [~] = optimize(cons,objs,opts);
    out.c = value(center_);
    out.G = out.G*diag(value(phi_));
end
end