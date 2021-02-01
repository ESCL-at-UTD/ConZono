function out = rescale(obj)

obj0 = copy(obj);

[~,E,~] = bounds(obj,1);

obj0.c = obj.c + obj.G*(E(:,2)+E(:,1))/2;
obj0.G = obj.G*diag(E(:,2)-E(:,1))/2;
obj0.b = obj.b - obj.A*(E(:,2)+E(:,1))/2;
obj0.A = obj.A*diag(E(:,2)-E(:,1))/2;

out = obj0;
end