function [out_R, out_E, out_rep] = bounds_ind(obj,maxIter)


e_lb = -ones(obj.nG,obj.nC); % Initialization
e_ub = ones(obj.nG,obj.nC);

E_old = [e_lb e_ub];
for out_rep = 1:maxIter
    r_lb = -inf*ones(obj.nG,obj.nC); % Initialization 
    r_ub =  inf*ones(obj.nG,obj.nC);

    for i = 1:obj.nC % Traverses through the constraints
        for j = 1:obj.nG % Traverses through the generators
            if obj.A(i,j) >= eps % Checks if A(i,j) \neq 0
                rj_temp = [0 0];
                for k = 1:obj.nG
                    if k ~= j
                        if obj.A(i,k)/obj.A(i,j) >= 0 % Computes bounds using interval arithmetic
                            rj_temp(1) = rj_temp(1) + obj.A(i,k)/obj.A(i,j)*e_lb(k,i);
                            rj_temp(2) = rj_temp(2) + obj.A(i,k)/obj.A(i,j)*e_ub(k,i);
                        else
                            rj_temp(1) = rj_temp(1) + obj.A(i,k)/obj.A(i,j)*e_ub(k,i);
                            rj_temp(2) = rj_temp(2) + obj.A(i,k)/obj.A(i,j)*e_lb(k,i);
                        end
                    end
                end
                try
                rj_temp = obj.b(i)/obj.A(i,j) - [rj_temp(2) rj_temp(1)];              
                catch
                   disp(i) 
                end
                r_lb(j,i) = max(r_lb(j,i),rj_temp(1)); % Computes the intersected interval between r_lb and computed bound.
                r_ub(j,i) = min(r_ub(j,i),rj_temp(2)); % Computes the intersected interval between r_ub and computed bound.
                e_lb(j,i) = max(e_lb(j,i),r_lb(j,i));  % Computes the intersected interval between e_lb and computed bound.
                e_ub(j,i) = min(e_ub(j,i),r_ub(j,i));  % Computes the intersected interval between e_ub and computed bound.
            end
        end
    end
    out_R = [r_lb r_ub];
    out_E = [e_lb e_ub];
    if (obj.nC ~= 0)
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

end