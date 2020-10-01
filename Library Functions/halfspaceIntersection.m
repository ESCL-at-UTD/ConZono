function [x] = halfspaceIntersection(z,Y)

% halfspaceIntersection computes the intersection of zonotope and halfspace Y.

% Inputs: z - Zonotope in CG-Rep (as a struct variable with c,
% G) parameters satisfying z = {c +G\xi, ||\xi||_{\infty} \leq 1}. 
% Y - Halfspace (Polyhedron object in Multi-parametric toolbox (MPT)).
% satisfying Y = \{y, Py \leq q \}

% Returns a constrained zonotope (x) in CG-Rep as a struct variable with c,G,A,b
% c-center of zonotope, G-generator matrix,A,b satisfying 
% x = {c +G\xi, ||\xi||_{\infty} <= 1, A\xi = b (constraints)}

x = z;

for j = 1:size(Y.H,1)
    if abs(Y.H(j,end)-Y.H(j,1:end-1)*x.c) < sum(abs(Y.H(j,1:end-1)*x.G)) % Checks for zonotope intersection

        d_max = Y.H(j,end)-Y.H(j,1:end-1)*x.c+sum(abs(Y.H(j,1:end-1)*x.G));
        
        x.A = [x.A zeros(size(x.A,1),1); Y.H(j,1:end-1)*x.G d_max/2]; % Lines 21-24 adds the constraints
        x.b = [x.b; Y.H(j,end)-Y.H(j,1:end-1)*x.c-d_max/2];
        x.c = x.c;
        x.G = [x.G, zeros(size(x.G,1),1)];
    end
end
end
