classdef conZono < matlab.mixin.Copyable
    % Constrained zonotope of the form
    % X = { c + G \xi | ||\xi||_inf <= 1, A \xi = b }
    
    properties
        c  % Center vector
        G  % Generator matrix
        A  % Constraint matrix
        b  % Constraint vector
        n  % Dimension
        nG % Number of generators
        nC % Number of constraints
    end
    
    methods
        function obj = conZono
            % Construct a constrained zonotope
        end
        
        % Methods in separate files
        getDimensions(obj)
        plot(obj,varargin) % Works
        out = volume(obj)
        out = mtimes(obj1,obj2) % Works
        out = plus(obj1,obj2)
        out = halfspaceIntersection(obj,H)
        out = generalizedIntersection(obj1,obj2,R)
        out = convexHull(obj1,obj2)
        out = conZono2AHPoly(obj)
        out = containCheck(obj1,obj2)
		out = zonContainCheck(obj1,obj2,cons_in)
        out = hausdorffDistance(obj1,obj2)
        out = innnerApprox(obj1,obj2,varargin)
        out = innerApproxPointContain(obj1,obj2,method,point)
		[out_R, out_E, out_rep] = bounds_ind(obj,maxIter)
        [out_R, out_E, out_rep] = bounds(obj,maxIter)
		out = redundancy_indices(obj,maxIter)	
		[out, out_ind] = rrefcp(obj)
		out = removerowicolumnj(obj,i,j)
		out = pontryagin_diff(obj1,obj2,varargin);
		out = pontryagin_diff_1step(obj1,obj2);
		out = bounds_lp(obj)
		out = rescale(obj)
		out = volumeratio(obj1,obj2)
    end
end

