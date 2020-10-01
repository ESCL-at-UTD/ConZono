function [x] = generalizedIntersection(z,y,R)
% generalizedIntersection computes the generalized intersection over a matrix between
% zonotopes/constrained zonotopes 

% Inputs: z,y - Zonotope/Constrained zonotope in CG-Rep (as a struct variable with c,
% G,A,b) parameters satisfying z = {c +G\xi, ||\xi||_{\infty} \leq 1, A\xi
% = b}. Cases like when z is a zonotope and y is a constrained zonotope or vice-versa also works.

% Generalized intersection (z \cap_R y) satisfies Rz \in y. 
% Setting R = eye(n) where n is the dimensional space computes the normal
% intersection.

% Returns a constrained zonotope (x) in CG-Rep as a struct variable with c,G,A,b
% c-center of zonotope, G-generator matrix,A,b satisfying 
% x = {c +G\xi, ||\xi||_{\infty} < 1, A\xi = b (constraints)}

x.c = z.c;
x.G = [z.G, zeros(size(z.G,1),size(y.G,2))]; % Lines 18-21 computes the generalized intersection
x.A = blkdiag(z.A,y.A);  
x.A = [x.A;[R*z.G -y.G]];
x.b = [z.b;y.b;y.c-R*z.c];

end

