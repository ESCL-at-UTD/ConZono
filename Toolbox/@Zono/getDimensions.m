function getDimensions(obj)
% Identify the dimension, number of generators, and number of
% constraints
obj.n  = size(obj.c,1);
obj.nG = size(obj.G,2);
obj.A = zeros(0,obj.nG);
end