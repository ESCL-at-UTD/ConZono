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
    
    needsProjection = true;
    
    % input checking
    if isempty(varargin)
        if obj.n > 3
            error(['Cannot plot sets in ',num2str(obj.n),' dimensions! Please, specify dimensions to project into.']);
        end
        color = 'r'; alpha = 1;
        needsProjection = false;
        
    elseif length(varargin) == 1
        dims = varargin{1};
        color = 'r'; alpha = 1;
        
    elseif length(varargin) == 2
        if obj.n > 3
            error(['Cannot plot sets in ',num2str(obj.n),' dimensions! Please, specify dimensions to project into.']);
        end
        
        color = varargin{1}; alpha = varargin{2};
        needsProjection = false;
        
    elseif length(varargin) == 3
        dims = varargin{1}; color = varargin{2}; alpha = varargin{3};
    end
    
    
    if obj.nC == 0
        Box = Polyhedron('lb',-ones(obj.nG,1),'ub',ones(obj.nG,1));
        if needsProjection
            %Projects center and generators on the selected dimensions
            c = obj.c(dims);
            G = obj.G(dims,:);
            needsProjection = false;
        else
            c = obj.c;
            G = obj.G;
        end
    else
        
        Box = Polyhedron('lb',-ones(obj.nG,1),'ub',ones(obj.nG,1),'He',[obj.A obj.b]);
        if needsProjection
            %Projects center and generators on the selected dimensions
            c = obj.c(dims);
            G = obj.G(dims,:);
            needsProjection = false;
        else
            c = obj.c;
            G = obj.G;
        end
    end
    
    P = c + G*Box;
    
    if needsProjection
        plot(P.projection(dims),'color',color,'alpha',alpha);
    else
        plot(P,'color',color,'alpha',alpha);
    end
    
end