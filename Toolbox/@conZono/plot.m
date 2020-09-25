function plot(obj,varargin)
% Plot a constrained zonotope
% Optional arguments:
% dims: vector with dimentions to project and plot the object
% color: color definition of the interior of the object to be plotted
% alpha: transparency of the color defined
% Default color: red. Default alpha: 1. Default dims: [1 2 3], [1 2], [1]
% for 3-, 2- and 1-dimensional objects respectively.
% 
% Usage options:
% obj.plot([1 2])
% obj.plot('b',0.1)
% obj.plot([1 2],'b',0.1);
%

obj.getDimensions;
% input checking
if isempty(varargin)
    if obj.n > 3
        error(['Cannot plot sets in ',num2str(obj.n),' dimensions!']);
    else
            color = 'r'; alpha = 1;
            if obj.n == 3
                dims = [1 2 3];
            elseif obj.n ==2
                dims = [1 2];
            else
                dims = 1;
            end
    end
elseif length(varargin) == 1
    dims = varargin{1};
    if ~(isinteger(dims) && size(dims,1)==1)
        error('Invalid input arguments.');
    end
    
elseif length(varargin) == 2
    color = varargin{1}; alpha = varargin{2};
      if obj.n == 3
                dims = [1 2 3];
            elseif obj.n ==2
                dims = [1 2];
            else
                dims = 1;
      end
elseif length(varargin) == 3
     dims = varargin{1}; color = varargin{2}; alpha = varargin{3}; 
      if ~(isinteger(dims) && size(dims,1)==1)
        error('Invalid input arguments.');
    end
end


    if obj.nC == 0
        Box = Polyhedron('lb',-ones(obj.nG,1),'ub',ones(obj.nG,1));
    else
        Box = Polyhedron('lb',-ones(obj.nG,1),'ub',ones(obj.nG,1),'He',[obj.A obj.b]);
    end
    P = plus(obj.c,affineMap(Box,obj.G));
    plot(P.projection(dims),'color',color,'alpha',alpha);

end