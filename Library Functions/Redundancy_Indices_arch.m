function [redund] = Redundancy_Indices(x,maxIter)
% Redundancy_Indices identifies indices of redundant constraints and generators of a constrained zonotope.

% Inputs - x (Constrained zonotope) in CG-Rep satisfies r = {c + G\xi, ||\xi||_{\infty} \leq 1, A\xi = b}
% maxIter - Maximum number of iterations to refine the interval set $E$ in function Bounds_ind

% Returns redund - matrix of redundant constraint and generator indices of identified redundancy

nc = size(x.A,1);
ng = size(x.A,2);
[R_ind,~,~] = Bounds_ind(x,maxIter); % Compute the Bounds of the generators for each constraint separately
redund = [];
for i = 1:nc
    for j = 1:ng
        if (R_ind(j,i) >= -1) && (R_ind(j,i+nc) <= 1)
            redund = [redund;i j];
        end
    end
end

end