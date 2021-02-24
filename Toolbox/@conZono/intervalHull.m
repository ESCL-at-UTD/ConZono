function [IH,IH_bounds] = intervalHull(obj)
% Computes the interval hull of a constrained zonotope
% Inputs:
%   - obj: conZono
% Outputs:
%   - IH: interval hull of 'obj' as a conZono
%   - IH_bounds: 2D matrix with lower and upper bounds that define the
%   interval hull (one row for each dimension)
% 
% Example:
% ng = 10;  
% nc = ng/2;
% c = zeros(10,1);
% G = rand(10,ng);
% A = rand(nc,ng);
% b = rand(nc,1);
% Z = conZono(c,G,A,b);
% [IH,IH_bounds] = Z.intervalHull;
%

%check for empty constrained zonotope
if obj.n==0 
    IH = conZono;
    IH_bounds = [];
    return;
end

%from eq. (3) of "Methods for Order Reduction of Zonotopes", Anna-Kathrin Kopetzki, Bastian Schuermann, and Matthias Althoff
% this is just to find the limits using the unconstrained zonotope, which
% ignores the A and b matrices of the constrained zonotope
% the interval [c-G_sum, c+G_sum] is the interval hull of the unsconstrained zonotope and
% thus overapproximates the constrained zonotope
G_sum = sum(abs(obj.G),2);
IH_lb = obj.c - G_sum;
IH_ub = obj.c + G_sum;


%unconstrained zonotope
if isempty(obj.A)
    IH_bounds=[IH_lb IH_ub];
else
 
    
    
    % x is the point being sought
    % \xi are the variables from the generators
    % The constraints are:
    % x = G\xi + c ==> x-G\xi = c ==> [I -G][x \xi]' = c
    % A\xi = b
    % Putting it all together: [I -G;0 A][x \xi]' = [c;b]
    n = obj.n; %num dimensions
    Ae = [eye(n) -obj.G;
        zeros(obj.nC,n) obj.A];
    be = [obj.c;obj.b];
    
    %bounds for the decision variables of the optimization problems below
    lb = [IH_lb;-ones(size(obj.A,2),1)];
    ub = [IH_ub;ones(size(obj.A,2),1)];
    
    IH_bounds = zeros(n,2);
    for i=1:n
        %we want to minimize f*[x \xi]
        f = zeros(1,n+obj.nG);
        f(i) = 1; %negative direction of dimension i
        lp = Opt('f',f,'Ae',Ae,'be',be,'lb',lb,'ub',ub); % Formulates the linear program
        sol = mpt_solve(lp); % Solves the LP.
        if sol.exitflag == 1
            IH_bounds(i,1) = sol.xopt(i);
        else
            error("Error finding solution to the LP. i=%d",i);
        end
        
        %we want to minimize f*[x \xi]
        f = zeros(1,n+obj.nG);
        f(i) = -1; %positive direction of dimension i
        lp = Opt('f',f,'Ae',Ae,'be',be,'lb',lb,'ub',ub); % Formulates the linear program
        sol = mpt_solve(lp); % Solves the LP.
        if sol.exitflag == 1
            IH_bounds(i,2) = sol.xopt(i);
        else
            error("Error finding solution to the LP. i=%d",i);
        end
    end
    
end

%creating the constrained zonotope based on the computed bounds
IH_lb = IH_bounds(:,1);
IH_ub = IH_bounds(:,2);
IH_c = IH_lb + (IH_ub-IH_lb)/2;
IH = conZono(IH_c,diag(IH_ub-IH_c));

end