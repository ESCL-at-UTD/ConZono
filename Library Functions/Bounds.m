function [R,E,rep] = Bounds(x,maxIter)

% Bounds computes the R, E bounds of the \xi variables associated with
% each of the generators of x considering all constraints.

% Inputs: x - constrained zonotope in CG-Rep as a struct variable with c,G,A,b
% c-center of zonotope, G-generator matrix,A,b satisfying 
% x = {c +G\xi, ||\xi||_{\infty} <= 1, A\xi = b (constraints)}
% maxIter - Maximum number of iterations for refinement of intervals R, E

% Returns intervals R_j \equiv [R_lb(j) R_ub(j)], E_j \equiv [E_lb(j) E_ub(j)]
% associated with each of the generators of x, and rep (number of iterations 
% needed for refinement of R, E)

nc = size(x.A,1);
ng = size(x.A,2);

e_lb = -ones(ng,1); % Initialization 
e_ub =  ones(ng,1);

E_old = [e_lb e_ub];

for rep = 1:maxIter
    r_lb = -inf*ones(ng,1); % Initialization 
    r_ub =  inf*ones(ng,1);

    for i = 1:nc % Traverses through the constraints
        for j = 1:ng % Traverses through the generators
            if x.A(i,j) >= eps % Checks if A(i,j) \neq 0
                rj_temp = [0 0];
                for k = 1:ng
                    if k ~= j
                        if x.A(i,k)/x.A(i,j) >= 0 % Computes bounds using interval arithmetic
                            rj_temp(1) = rj_temp(1) + x.A(i,k)/x.A(i,j)*e_lb(k);
                            rj_temp(2) = rj_temp(2) + x.A(i,k)/x.A(i,j)*e_ub(k);
                        else
                            rj_temp(1) = rj_temp(1) + x.A(i,k)/x.A(i,j)*e_ub(k);
                            rj_temp(2) = rj_temp(2) + x.A(i,k)/x.A(i,j)*e_lb(k);
                        end
                    end
                end                        
                rj_temp = x.b(i)/x.A(i,j) - [rj_temp(2) rj_temp(1)];
                r_lb(j) = max(r_lb(j),rj_temp(1)); % Computes the intersected interval between r_lb and computed bound among all the constraints.
                r_ub(j) = min(r_ub(j),rj_temp(2)); % Computes the intersected interval between r_ub and computed bound among all the constraints.
                e_lb(j) = max(e_lb(j),r_lb(j));    % Computes the intersected interval between e_lb and computed bound among all the constraints.
                e_ub(j) = min(e_ub(j),r_ub(j));    % Computes the intersected interval between e_ub and computed bound among all the constraints.
            end
        end
    end        
    R = [r_lb r_ub];
    E = [e_lb e_ub];
    sum_feas = 0;
    for j = 1:ng
        if (e_lb(j) <= e_ub(j)) % Checks if interval set E is feasible.
        sum_feas = sum_feas + 1;
        end
    end
    if sum_feas ~= ng
        break
    end
    if isequal(E, E_old) % Checks if the interval set E has been refined since the last iteration       
        break        
    else
        E_old = E;
    end
end
end
