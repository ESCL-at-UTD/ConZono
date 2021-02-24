function out = convexHull(obj1,obj2)

out = conZono;
if obj1.n ~= obj2.n
    disp(['Cannot compute the convex hull of a zonotope in ',num2str(obj1.n),...
        ' dimensions to a zonotope in ',num2str(obj2.n),' dimensions!'])
else
    out.c = (obj1.c + obj2.c)/2;
    out.G = [obj1.G obj2.G (obj1.c-obj2.c)/2 zeros(obj1.n,2*(obj1.nG+obj2.nG))];
    H = [eye(obj1.nG) zeros(obj1.nG,obj2.nG) -0.5*ones(obj1.nG,1);...
        -eye(obj1.nG) zeros(obj1.nG,obj2.nG) -0.5*ones(obj1.nG,1);...
        zeros(obj2.nG,obj1.nG) eye(obj2.nG)  0.5*ones(obj2.nG,1);...
        zeros(obj2.nG,obj1.nG) -eye(obj2.nG)  0.5*ones(obj2.nG,1)];
    I = 1*eye(2*(obj1.nG + obj2.nG));
    f = -0.5*ones(2*(obj1.nG + obj2.nG),1);
    if isempty(obj1.b)
        b1 = zeros(obj1.nC,1);
    else
        b1 = obj1.b;
    end
    if isempty(obj2.b)
        b2 = zeros(obj2.nC,1);
    else
        b2 = obj2.b;
    end
    out.A = [obj1.A zeros(obj1.nC,obj2.nG) -b1/2 zeros(obj1.nC,size(out.G,2)-obj1.nG-obj2.nG-1);...
        zeros(obj2.nC,obj1.nG) obj2.A +b2/2 zeros(obj2.nC,size(out.G,2)-obj2.nG -obj1.nG-1);...
        H I];
    out.b = [obj1.b/2;obj2.b/2;f];
end
end