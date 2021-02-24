function out = generalizedIntersection(obj1,obj2,R)

out = conZono;
if (obj1.n ~= size(R,2)) || (obj2.n ~= size(R,1))
    disp('Inconsistent dimensions!')
else
    out.c = obj1.c;
    out.G = [obj1.G, zeros(size(obj1.G,1),size(obj2.G,2))];
    out.A = blkdiag(obj1.A,obj2.A);
    out.A = [out.A;[R*obj1.G -obj2.G]];
    out.b = [obj1.b;obj2.b;obj2.c-R*obj1.c];
end
end