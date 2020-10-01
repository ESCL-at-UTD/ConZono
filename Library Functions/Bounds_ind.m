function [R,E,rep] = Bounds_ind(x,maxIter)
% Bounds_ind computes the R, E bounds of the \xi variables associated with
% each of the generators of x with each constraint taken seperately.
% Inputs: x - constrained zonotope in CG-Rep as a struct variable with c,G,A,b
% c-center of zonotope, G-generator matrix,A,b satisfying 
% x = {c +G\xi, ||\xi||_{\infty} <= 1, A\xi = b (constraints)}
% maxIter - Maximum number of iterations for refinement of intervals R, E

% Returns R_j \equiv [R_lb(j,i) R_ub(j,i)] and E_j \equiv [E_lb(j,i) E_ub(j,i)]
% associated with each of the generators of x, and rep (number of iterations 
% needed for refinement of R, E)

nc = size(x.A,1);
ng = size(x.A,2);

e_lb = -ones(ng,nc); % Initialization 
e_ub =  ones(ng,nc);

E_old = [e_lb e_ub];

for rep = 1:maxIter
    r_lb = -inf*ones(ng,nc); % Initialization 
    r_ub =  inf*ones(ng,nc);

    for i = 1:nc % Traverses through the constraints
        for j = 1:ng % Traverses through the generators
            if x.A(i,j) >= eps % Checks if A(i,j) \neq 0
                rj_temp = [0 0];
                for k = 1:ng
                    if k ~= j
                        if x.A(i,k)/x.A(i,j) >= 0 % Computes bounds using interval arithmetic
                            rj_temp(1) = rj_temp(1) + x.A(i,k)/x.A(i,j)*e_lb(k,i);
                            rj_temp(2) = rj_temp(2) + x.A(i,k)/x.A(i,j)*e_ub(k,i);
                        else
                            rj_temp(1) = rj_temp(1) + x.A(i,k)/x.A(i,j)*e_ub(k,i);
                            rj_temp(2) = rj_temp(2) + x.A(i,k)/x.A(i,j)*e_lb(k,i);
                        end
                    end
                end                        
                rj_temp = x.b(i)/x.A(i,j) - [rj_temp(2) rj_temp(1)];
                r_lb(j,i) = max(r_lb(j,i),rj_temp(1)); % Computes the intersected interval between r_lb and computed bound.
                r_ub(j,i) = min(r_ub(j,i),rj_temp(2)); % Computes the intersected interval between r_ub and computed bound.
                e_lb(j,i) = max(e_lb(j,i),r_lb(j,i));  % Computes the intersected interval between e_lb and computed bound.
                e_ub(j,i) = min(e_ub(j,i),r_ub(j,i));  % Computes the intersected interval between e_ub and computed bound.
            end
        end
    end
    R = [r_lb r_ub];
    E = [e_lb e_ub];
    if (nc ~= 0)
	sum_feas = 0; % Lines 50-58 additions
        for j = 1:ng
            if (e_lb(j) <= e_ub(j)) % Checks if interval set E is feasible.
            sum_feas = sum_feas + 1;
            end
        end
        if sum_feas ~= ng
            break
        end
    end
    if isequal(E, E_old)
        break
    else
        E_old = E;
    end
end
end

