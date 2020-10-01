function out = mtimes(obj,a)
obj.getDimensions;
if ~isscalar(a)
    disp('Can only multiply AH-polytopes by scalars!')
else
    out = affineHPoly;
    out.c = a*obj.c;
    out.G = obj.G;
    out.A = obj.A;
    out.b = a*obj.b;
    out.getDimensions;
end
end