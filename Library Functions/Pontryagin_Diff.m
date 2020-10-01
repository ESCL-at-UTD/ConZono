function [Zint] = Pontryagin_Diff(Z1,Z2,removeRedund,maxIter)
% Pontryagin_Diff computes the Pontryagin difference between zonotopes.
% Also works if Z_1 is a constrained zonotope.

% Inputs: Z1 - Zonotope/Constrained Zonotope in CG-Rep (as a struct variable with c_1,
% G_1, A_1, b_1) parameters satisfying Z1 = {c_1 +G_1\xi_1,
% ||\xi_1||_{\infty} \leq 1, A_1\xi_1 = b_1}. 
% Z2 - Zonotope in CG-Rep (as a struct variable with c_2, G_2) parameters 
% satisfying Z2 = {c_2 +G_2\xi_2, ||\xi_2||_{\infty} \leq 1}. 

% Returns a constrained zonotope (Zint) in CG-Rep as a struct variable with
% c_d, G_d, A_d, b_d
% c_d, G_d, A_d, b_d satisfying 
% x_d = {c_d +G_d\xi_d, ||\xi_d||_{\infty} < 1, A_d\xi_d = b_d (constraints)}

n_dof = size(Z2.G,1);
ng_sub = size(Z2.G,2);

Zint = [];

Zint.c = Z1.c - Z2.c;
Zint.G = Z1.G;
Zint.A = Z1.A;
Zint.b = Z1.b;

for i =1 : ng_sub % Iterates over the generators of $Z_2$
    Zint_plus.c = Zint.c + Z2.G(:,i); % Shifts Zint by Z2 to the right
    Zint_plus.G = [Zint.G];
    Zint_plus.A = [Zint.A]; 
    Zint_plus.b = Zint.b;
    
    Zint_minus.c = Zint.c - Z2.G(:,i); % Shifts Zint by Z2 to the left
    Zint_minus.G = Zint.G ;
    Zint_minus.A = Zint.A;                 
    Zint_minus.b = Zint.b;
    
    [Zint] = generalizedIntersection(Zint_plus,Zint_minus,eye(n_dof)); % Shifted zonotope intersection
    if removeRedund == 1 % Remove redundancy
        [redund] = Redundancy_Indices(Zint,maxIter);
        while size(redund,1) > 0
            [Zint] = RemoveRowiColumnj(Zint,redund(1,1),redund(1,2));
            [redund] = Redundancy_Indices(Zint,maxIter);
        end
    end
end

end