function R = reachDT(A,B,X0,U,N,Xfeasible)
    % Computes the forward reachable sets for a discrete time linear system of
    % the form x+ = Ax + Bu;
    %
    % Inputs:
    % A,B: discrete time system matrices
    % X0: initial set
    % U: input set
    % N: number of steps to compute the reachable sets
    % Xfeasible: a Polyhedron. If passed, each reachable set is intersected with the set Xfeasible
    % before computing the next reachable set. X0 is also intersected with
    % Xfeasible. In this case, X0 and U must be conZono class instances.
    %
    % Outputs:
    % R: cell array of reachable sets
    %
    % Example:
    % Ac = [-1 -4; 4 -1];
    % Bc = [1;0];
    % sys_c = ss(Ac,Bc,eye(2),zeros(2,1));
    % dt = 0.1;
    % sys_d = c2d(sys_c,dt);
    % A = sys_d.A;
    % B = sys_d.B;
    %
    % % Initial set for reachability
    % X0 = conZono;
    % X0.c = [5;5];
    % X0.G = eye(2);
    %
    % % Input set
    % U = conZono;
    % U.c = 0;
    % U.G = 1;
    %
    % N=10;
    % R = reachDT(A,B,X0,U,N);
    
    
    doIntersection = false;
    
    %input arguments checks
    if(~isa(X0,'conZono') && ~isa(X0,'Zono') && ~isa(X0,'affineHPoly'))
        error('Argument X0 must be a conZono, Zono or affineHPoly');
    elseif (~isa(U,'conZono') && ~isa(U,'Zono') && ~isa(U,'affineHPoly'))
        error('Argument U must be a conZono, Zono or affineHPoly');
    elseif (nargin==6)
        doIntersection = true;
        if (~isa(Xfeasible,'Polyhedron'))
            error('Argument Xfeasible must be a Polyhedron');
        elseif (~isa(X0,'conZono') && ~isa(U,'Zono'))
            error('If argument Xfeasible is provided, X0 and U must be conZono class instances');
        end
    end
    
    R = cell(N+1,1);
    
    R{1} = X0;
    if doIntersection
        R{1} = halfspaceIntersection(R{1},Xfeasible);
    end
    
    BU = B*U;
    %compute reachable sets
    for i=1:N
        R{i+1} = A*R{i} + BU;
        
        %intersects with XFeasible, if supplied
        if doIntersection
            R{i+1} = halfspaceIntersection(R{i+1},Xfeasible);
        end
    end
    
    
    
    
    
    