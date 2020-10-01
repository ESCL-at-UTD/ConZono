function [y] = zono_minHRep(x)
% zono_minHRep converts the zonotope to minimal H-Rep 

% Inputs: x (zonotope) (in CG-Rep) (struct variable with c,G parameters) 
% c-center of zonotope, G-generator matrix,
% x = {c +G\xi, ||\xi||_{\infty} \leq 1}

% Returns a Polyhedron in minimal H-Rep.
    r = x;

    Box_r = Polyhedron('lb',-ones(size(r.G,2),1),'ub',ones(size(r.G,2),1));
    y = plus(affineMap(Box_r,r.G),r.c);
end