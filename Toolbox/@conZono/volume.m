function out = volume(obj)
if obj.nC == 0
    Box = Polyhedron('lb',-ones(obj.nG,1),'ub',ones(obj.nG,1));
else
    Box = Polyhedron('lb',-ones(obj.nG,1),'ub',ones(obj.nG,1),'He',[obj.A obj.b]);
end
P = plus(obj.c,affineMap(Box,obj.G));
out = P.volume;
end