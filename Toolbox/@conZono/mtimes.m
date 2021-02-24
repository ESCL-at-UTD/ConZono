function out = mtimes(obj1,obj2)
if strcmp(class(obj1),'conZono')
    obj = copy(obj1);
    a = obj2;
elseif strcmp(class(obj2),'conZono')
    obj = copy(obj2);
    a = obj1;
else
    disp('Can only multiply zonotopes by scalars or matrices!')
    return
end

out = conZono;
out.c = a*obj.c;
out.G = a*obj.G;
out.A = obj.A;
out.b = obj.b;

end