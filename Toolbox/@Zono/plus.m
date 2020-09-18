function out = plus(obj1,obj2)
obj1.getDimensions;
obj2.getDimensions;
if obj1.n ~= obj2.n
    disp(['Cannot add a zonotope in ',num2str(obj1.n),...
        ' dimensions to a zonotope in ',num2str(obj2.n),' dimensions!'])
else
    out = conZono;
    out.c = obj1.c + obj2.c;
    out.G = [obj1.G obj2.G];
    out.getDimensions;
end
end