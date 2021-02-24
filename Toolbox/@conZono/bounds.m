function [out_R,out_E,out_rep] = bounds(obj,maxIter)


e_lb = -ones(obj.nG,1); % Initialization 
e_ub =  ones(obj.nG,1);

E_old = [e_lb e_ub];

for out_rep = 1:maxIter
    r_lb = -inf*ones(obj.nG,1); % Initialization 
    r_ub =  inf*ones(obj.nG,1);

    for i = 1:obj.nC % Traverses through the constraints
        for j = 1:obj.nG % Traverses through the generators
            if obj.A(i,j) >= eps % Checks if A(i,j) \neq 0
                rj_temp = [0 0];
                for k = 1:obj.nG
                    if k ~= j
                        if obj.A(i,k)/obj.A(i,j) >= 0 % Computes bounds using interval arithmetic
                            rj_temp(1) = rj_temp(1) + obj.A(i,k)/obj.A(i,j)*e_lb(k);
                            rj_temp(2) = rj_temp(2) + obj.A(i,k)/obj.A(i,j)*e_ub(k);
                        else
                            rj_temp(1) = rj_temp(1) + obj.A(i,k)/obj.A(i,j)*e_ub(k);
                            rj_temp(2) = rj_temp(2) + obj.A(i,k)/obj.A(i,j)*e_lb(k);
                        end
                    end
                end                        
                rj_temp = obj.b(i)/obj.A(i,j) - [rj_temp(2) rj_temp(1)];
                r_lb(j) = max(r_lb(j),rj_temp(1)); % Computes the intersected interval between r_lb and computed bound among all the constraints.
                r_ub(j) = min(r_ub(j),rj_temp(2)); % Computes the intersected interval between r_ub and computed bound among all the constraints.
                e_lb(j) = max(e_lb(j),r_lb(j));    % Computes the intersected interval between e_lb and computed bound among all the constraints.
                e_ub(j) = min(e_ub(j),r_ub(j));    % Computes the intersected interval between e_ub and computed bound among all the constraints.
            end
        end
    end        
    out_R = [r_lb r_ub];
    out_E = [e_lb e_ub];
    sum_feas = 0;
    for j = 1:obj.nG
        if (e_lb(j) <= e_ub(j))
        sum_feas = sum_feas + 1;
        end
    end
    if sum_feas ~= obj.nG
        break
    end
    if isequal(out_E, E_old)        
        break        
    else
        E_old = out_E;
    end
end

end