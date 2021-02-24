function r = Rescale(x,E)

% Rescale generates the exact set (x) but with different bounds of xi variables extending from -1 to 1.

% Inputs: x - constrained zonotope in CG-Rep as a struct variable with c,G,A,b
% c-center of zonotope, G-generator matrix,A,b satisfying 
% x = {c +G\xi, ||\xi||_{\infty} \leq 1, A\xi = b (constraints)}
% E - Interval set capturing the bounds of xi variables of all the generators

% Returns re-scaled constrained zonotope r in CG-Rep as a struct variable with c_r,G_r,A_r,b_r
% c_r-center of zonotope, G_r-generator matrix,A_r,b_r satisfying 
% x = {c_r +G_r\xi_r, ||\xi_r||_{\infty} \leq 1, A_r\xi_r = b_r (constraints)}

r.c = x.c + x.G*(E(:,2)+E(:,1))/2;
r.G = x.G*diag(E(:,2)-E(:,1))/2;
r.b = x.b - x.A*(E(:,2)+E(:,1))/2;
r.A = x.A*diag(E(:,2)-E(:,1))/2;

end