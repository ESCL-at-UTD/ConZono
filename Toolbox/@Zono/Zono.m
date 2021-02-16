classdef Zono < matlab.mixin.Copyable
    % Zonotope of the form
    % X = { c + G \xi | ||\xi||_inf <= 1 }
    
    properties
        c  % Center vector
        G  % Generator matrix 
    end
    
    %These properties values depend on other properties and do not store any data themselves
    properties (Dependent)
        n  % Dimension
        nG % Number of generators
        order % Zonotope order
    end
    
    methods
        function obj = Zono
            % Construct a zonotope
        end
        
         %Property getter methods
        function value = get.n(obj)
            value  = size(obj.c,1);
        end
        function value = get.nG(obj)
            value  = size(obj.G,2);
        end
        function value = get.order(obj)
            value  = obj.nG/obj.n;
        end
        
        % Methods in separate files
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

