classdef conZono < matlab.mixin.Copyable
    % Constrained zonotope of the form
    % X = { c + G \xi | ||\xi||_inf <= 1, A \xi = b }
    
    properties
        c  % Center vector
        G  % Generator matrix
        A  % Constraint matrix
        b  % Constraint vector
    end
    
    %These properties values depend on other properties and do not store any data themselves
    properties (Dependent)
        n  % Dimension
        nG % Number of generators
        nC % Number of constraints
        order % Constrained zonotope order
    end
    
    methods
        function obj = conZono(c,G,A,b)
            if nargin==0
                %empty constructor
            elseif nargin==2
                % Construct a zonotope as a constrained zonotope providing
                % c and G only
                obj.c = c;
                obj.G = G;
            elseif nargin == 4
                 % Construct a constrained zonotope providing all matrices
                obj.c = c;
                obj.G = G;
                obj.A = A;
                obj.b = b;            
            else
                error('Invalid number of arguments for conZono constructor.');
            end
        end
        
        
         function value = get.A(obj)
            
            %This code ensures that, if A is empty, any code accessing this property will
            %get an empty 0 by obj.nG matrix. This ensures consistency in several
            %methods when concatenating A with other matrices without
            %having to check if A is empty.
            if size(obj.A,1)==0
                value = zeros(0,obj.nG);
            else
                value = obj.A;
            end
            
        end
        
        
        %Property getter methods
        function value = get.n(obj)
            value  = size(obj.c,1);
        end
        function value = get.nG(obj)
            value  = size(obj.G,2);
        end
        function value = get.nC(obj)
            value  = size(obj.A,1);
        end
        function value = get.order(obj)
            value  = (obj.nG - obj.nC)/obj.n;
        end
        
        
        % Methods in separate files
        plot(obj,varargin) 
        out = volume(obj,use_old_method)
        out = mtimes(obj1,obj2)
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
        out = cartProd(obj1,obj2)
        [IH,IH_bounds] = intervalHull(obj)
    end
end

