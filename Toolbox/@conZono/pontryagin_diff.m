function [out] = pontryagin_diff(obj1,obj2,varargin)

if obj1.n ~= obj2.n
	disp(['Cannot compute the Pontryagin difference between a zonotope in ',num2str(obj1.n),...
        ' dimensions and a zonotope in ',num2str(obj2.n),' dimensions!'])
elseif obj2.nC ~= 0
	disp(['Cannot compute the Pontryagin difference between a zonotope in ',num2str(obj1.n),...
        ' dimensions and a constrained zonotope'])
else 

    if isempty(varargin)
        removeRedund = 1;
		maxIter = 100;
    else
        removeRedund = varargin(1);
		maxIter = varargin(2);
    end

out = conZono;

out.c = obj1.c - obj2.c;
out.G = obj1.G;
out.A = obj1.A;
out.b = obj1.b;
out_plus = conZono;
out_minus = conZono;

for i =1 : obj2.nG
    out_plus.c = out.c + obj2.G(:,i);
    out_plus.G = [out.G];
    out_plus.A = [out.A]; 
    out_plus.b = out.b;
    
    out_minus.c = out.c - obj2.G(:,i);
    out_minus.G = out.G;
    out_minus.A = out.A;                 
    out_minus.b = out.b;
    
    [out] = generalizedIntersection(out_plus,out_minus,eye(obj1.n));
    if removeRedund == 1
        [redund] = redundancy_indices(out,maxIter);
        while size(redund,1) > 0
            [out] = removerowicolumnj(out,redund(1,1),redund(1,2)); % To add 
            [redund] = redundancy_indices(out,maxIter);
        end
    end
end

end

end