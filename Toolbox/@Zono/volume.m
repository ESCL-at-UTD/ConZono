function out = volume(obj)
obj.getDimensions;
Box = Polyhedron('lb',-ones(obj.nG,1),'ub',ones(obj.nG,1));
P = plus(obj.c,affineMap(Box,obj.G));
out = P.volume;
end