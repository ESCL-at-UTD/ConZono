classdef Zono < matlab.mixin.Copyable
    % Zonotope of the form
    % X = { c + G \xi | ||\xi||_inf <= 1 }
    
    properties
        c  % Center vector
        G  % Generator matrix
        A  % Constraint matrix
        b  % Constraint vector
        n  % Dimension
        nG % Number of generators   
    end
    
    methods
        function obj = Zono
            % Construct a zonotope
        end
        
        % Methods in separate files
        getDimensions(obj)
        plot(obj,varargin)
		out = containCheck(obj1,obj2)
        out = innnerApprox(obj1,obj2,varargin)
		out = redundancy_parallelgen(obj,varargin)		
		out = Zono2AHPoly(obj)
		out = volume(obj)
		out = volumeratio(obj1,obj2)		
		out = plus(obj1,obj2)
    end
end

