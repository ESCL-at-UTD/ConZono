function out = cartProd(obj1,obj2)
% Computes the cartesian product of two constrained zonotopes
% Inputs:
%   - obj1: conZono
%   - obj2: conZono
% Outputs:
%   - out: cartesian product (obj1 X obj2) as a conZono
% 
% Example:
% z1 = conZono([10;10],[1;1],[3 4; -5 2],[0;0]);
% z2 = conZono([0],[5]);
% 
% z3 = cartProd(z1,z2);
%
% -----------------------------------------------------------------------------------------------------------------
% 2021-01-12: Adapted from CORA 2020 (https://github.com/TUMcps/CORA) 'cartProd' function for constrained zonotopes
% -----------------------------------------------------------------------------------------------------------------
%

arguments
   obj1
   obj2 conZono
end


% new center vector
c = [obj1.c; obj2.c];

% new generator matrix
G = blkdiag(obj1.G,obj2.G);

% new constraint matrix
if isempty(obj1.A)
    if isempty(obj2.A)
        A = [];
    else
        A = [zeros(obj2.nC,obj1.nG),obj2.A];
    end
else
    if isempty(obj2.A)
        A = [obj1.A,zeros(obj1.nC,obj2.nG)];
    else
        A = [[obj1.A,zeros(obj1.nC,obj2.nG)];[zeros(obj2.nC,obj1.nG),obj2.A]];
    end
end

% new constraint offset
b = [obj1.b;obj2.b];

% generate resulting constrained zonotope
out = conZono(c,G,A,b);

end

