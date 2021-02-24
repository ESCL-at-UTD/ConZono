function out = volume(obj)
Box = Polyhedron('lb',-ones(obj.nG,1),'ub',ones(obj.nG,1));
P = plus(obj.c,affineMap(Box,obj.G));
out = P.volume;
end