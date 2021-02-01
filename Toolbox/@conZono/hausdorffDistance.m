function out = hausdorffDistance(obj1,obj2)

if obj1.n ~= obj2.n
    disp(['Cannot compute the Hausdorff distance of a zonotope in ',num2str(obj1.n),...
        ' dimensions to a zonotope in ',num2str(obj2.n),' dimensions!'])
else
    h1 = conZono2AHPoly(obj1);
    h2 = conZono2AHPoly(obj2);
    out = hausdorffDistance(h1,h2);
end
end