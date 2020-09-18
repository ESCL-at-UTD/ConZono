classdef affineHPoly < matlab.mixin.Copyable
    % Affine H-polytope (AH-polytope) of the form
    % X = { c + G \xi | A \xi <= b }
    
    properties
        c  % Center vector
        G  % Generator matrix
        A  % Constraint matrix
        b  % Constraint vector
        n  % Dimension
        nG % Number of generators
        nH % Number of halfspaces
    end
    
    methods
        function obj = affineHPoly
            % Construct a AH-polytope
        end
        
        % Methods in separate files
        getDimensions(obj)
        plot(obj,varargin)
        out = mtimes(obj,a)
        out = plus(obj1,obj2)
        out = containCheck(obj1,obj2)
        cons_out = conContainCheck(obj1,obj2,cons_in)
        out = hausdorffDistance(obj1,obj2)
		out = conPointContain(point,h0,cons_in);
    end
end

