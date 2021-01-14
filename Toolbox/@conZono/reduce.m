function reduced_obj = reduce(obj)

%TODO
% - include comments
% - include in examples

obj.getDimensions;

bound_iter = 100;

[R,~,~] = bounds(obj,bound_iter);
R_max = max(abs(R),[],2);
[~,j_min] = min(R_max);
[~,i_min] = max(abs(obj.A(:,j_min)));

reduced_obj = removerowicolumnj(obj,i_min,j_min);

reduced_obj.getDimensions;
end