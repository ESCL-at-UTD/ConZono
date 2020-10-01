function [ratio] = VolumeRatio(X,Y)
% VolumeRatio computes the volume ratio of two polytopes.

% Inputs: X - Polyhedron object in H-Rep
%         Y - Polyhedron object in H-Rep.

% Returns ratio = (X.volume/Y.volume)^(1/n) where $n$ is the dimension of X
% and Y.

n = size(X.H,2)-1;

ratio = (X.volume/Y.volume)^(1/n);
end

