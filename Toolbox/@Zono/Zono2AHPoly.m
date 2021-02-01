function out = Zono2AHPoly(obj)
out = affineHPoly;
out.c = obj.c;
out.G = obj.G;
out.A = [-eye(obj.nG);eye(obj.nG)];
out.b = ones(2*obj.nG,1);
end