function out = containCheck(obj1,obj2)
% Check if obj1 is contained in obj2

if obj1.n ~= obj2.n
    disp(['Cannot check containment of a zonotope in ',num2str(obj1.n),...
        ' dimensions to a zonotope in ',num2str(obj2.n),' dimensions!'])
else
    h1 = conZono2AHPoly(obj1);
    h2 = conZono2AHPoly(obj2);
    out = containCheck(h1,h2);
end
end