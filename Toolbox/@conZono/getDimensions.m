function getDimensions(obj)
% Identify the dimension, number of generators, and number of
% constraints
obj.n  = size(obj.c,1);
obj.nG = size(obj.G,2);
obj.nC = size(obj.A,1);
if obj.nC == 0
    obj.A = zeros(0,obj.nG);
end
end