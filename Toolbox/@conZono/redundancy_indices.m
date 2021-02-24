function [out] = redundancy_indices(obj,maxIter)

[R_ind,~,~] = bounds_ind(obj,maxIter);
out = [];

for i = 1:obj.nC
    for j = 1:obj.nG
        if (R_ind(j,i) >= -1) && (R_ind(j,i+obj.nC) <= 1)
            out = [out;i j];
        end
    end
end

end