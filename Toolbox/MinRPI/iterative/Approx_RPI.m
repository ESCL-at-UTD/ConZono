function [F_approx,s,alpha_o] = Approx_RPI(A,W,epsilon,MAX_ITER)
    arguments
        A
        W
        epsilon  (1,1)  double
        MAX_ITER (1,1)  {mustBeInteger}   = 100
        
    end
    
    if ~isprop(W,'H') && ~isfield(W,'G') %&& ~isprop(W,'G')
        error('Set must be in either H-Rep or G-Rep');
    end
    
    %%
    
    % Choose initial s
    s = 0;
    % Set maximum iterations
    for i = 1:MAX_ITER
        % Increment s
        s = s + 1;
        % Compute alpha min of s
        alpha_o = alpha_min(W,A,s);
        % Compute M(s)
        M_s = M_of_s(W,A,s,MAX_ITER);
        if alpha_o <= epsilon/(epsilon+M_s)
            break
        end
    end
    
    Fs = W;
    
    
    
    %%
    if isprop(W,'H')
        
        for i = 1:s % because nilpotent in one step
            Fs = plus(Fs,affineMap(W,A^i));
            Fs.minHRep;
        end
        F_approx = 1/(1-alpha_o)*Fs;
        F_approx.minHRep;
        % F_approx.plot
    elseif isfield(W,'G') %|| isprop(W,'G')
        
        Ai = eye(size(A,1));
        for i = 1:s % because nilpotent in one step
            Ai = A*Ai; %to avoid recomputation of A^i at every step
            Fs.c = Fs.c + Ai*W.c;
            Fs.G = [Fs.G Ai*W.G]; %TO DO, this changes size at every iteration. Allocate before beginning of the for loop
        end
        F_approx.c = 1/(1-alpha_o)*Fs.c;
        F_approx.G = 1/(1-alpha_o)*Fs.G;
        
        F_approx = conZono(F_approx.c, F_approx.G); %converts to conZono
    else
        
    end
    
    