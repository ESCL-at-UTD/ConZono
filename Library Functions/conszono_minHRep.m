function [y] = conszono_minHRep(x,meth)
% conszono_minHRep converts the constrained zonotope to minimal H-Rep 

% Inputs: x (constrained zonotope) (in CG-Rep) (struct variable with c,G,A,b parameters) 
% c-center of zonotope, G-generator matrix,A,b satisfying 
% x = {c +G\xi, ||\xi||_{\infty} < 1, A\xi = b (constraints)}
% meth = 1/2 (Explicit Minkowski addition/ Variable substitution and Minkowski addition)

% Returns a Polyhedron in minimal H-Rep.
r = x;

if (meth == 1) % Computes minkowski sum explicitly by definition
    Box_r = Polyhedron('lb',-ones(size(r.G,2),1),'ub',ones(size(r.G,2),1),...
                 'He',[r.A r.b]);
    Box_r= minHRep(Box_r);
    y = plus(r.c,affineMap(Box_r,r.G));
    y = minHRep(y);
elseif (meth == 2) % Computes the minkowski sum implicitly using a variable substitution for \xi
    % Formulates \xi = T\xi^' + s (affine transformation)
    T = null(x.A,'r'); % Computes the right null space of A.
    s = pinv(x.A)*x.b; % Computes the pesudo-inverse
    P = Polyhedron('H',[T ones(size(s,1),1)-s; -T ones(size(s,1),1)+s]);
    P = minHRep(P); % Removes redundant halfspaces
    x_poly = plus(x.c+x.G*s,affineMap(P,x.G*T)); % Computes the constrained zonotope x in H-Rep
    y = minHRep(x_poly); % Removes redundant halfspaces
end
end