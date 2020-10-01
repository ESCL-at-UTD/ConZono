function [E] = Bounds_lp(x)

% Bounds_lp computes the E bounds of the \xi variables associated with
% each of the generators of x considering all constraints.
% Inputs: x - constrained zonotope in CG-Rep as a struct variable with c,G,A,b
% c-center of zonotope, G-generator matrix,A,b satisfying 
% x = {c +G\xi, ||\xi||_{\infty} < 1, A\xi = b (constraints)}

% Returns E_j \equiv [E_lb(j) E_ub(j)] associated with each of the generators of x.

ng = size(x.A,2);

for i = 1:ng
    s = zeros(ng,1);
    s(i) = 1;
    lp = Opt('f',s','Ae',x.A,'be',x.b,'lb',-ones(size(x.A,2),1),'ub',ones(size(x.A,2),1)); % Formulates the linear program
    opt = mpt_solve(lp); % Solves the LP.
    E(i,1) = s'*opt.xopt;
end

for i = 1:ng
    s = zeros(ng,1);
    s(i) = 1;
    lp = Opt('f',-s','Ae',x.A,'be',x.b,'lb',-ones(size(x.A,2),1),'ub',ones(size(x.A,2),1)); % Formulates the linear program
    opt = mpt_solve(lp); % Solves the LP.
    E(i,2) = s'*opt.xopt;
end

end

