function getDimensions(obj)
% Identify the dimension, number of generators, and number of halfspaces
obj.n  = size(obj.c,1);
obj.nG = size(obj.G,2);
obj.nH = size(obj.A,1);
end