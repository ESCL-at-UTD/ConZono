function [out] = redundancy_parallelgen(obj,varargin)

obj0 = copy(obj)

if isempty(varargin)
	eps = 1e-6;
end

ind_rem = [];
for i = 1: obj0.nG
	for j = i+1:obj0.nG
		bool_pl = abs(dot(obj0.G(:,i),obj0.G(:,j)))/(norm(obj0.G(:,i))*norm(obj0.G(:,j))) >= 1 - eps;
		if bool_pl == 1
            ind_rem = [ind_rem j]
			obj0.G(:,i) = obj0.G(:,i) + obj0.G(:,j);			
		end		
	end
end
obj0.G(:,ind_rem) = [];
out = obj0;
end