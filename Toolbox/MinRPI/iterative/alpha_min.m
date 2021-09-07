function alpha_o = alpha_min(W,A,s)

if isprop(W,'H')
    i_cons = size(W.H,1);
    
    alpha_i = zeros(i_cons,1);
    for i = 1:i_cons
        fi = W.H(i,1:end-1)';
        gi = W.H(i,end);
        a = (A^s)'*fi;
        hW_a = support(W,a);
        alpha_i(i) = hW_a/gi;
    end
    alpha_o = max(alpha_i);
elseif isfield(W,'G') || isprop(W,'G')
    i_gens = size(W.G,2);
    
    alpha_i = zeros(i_gens,1);
    for i = 1:i_gens
        gi = W.G(:,i);
        alpha_i(i) = ( sum(abs(gi'*A^s*W.G)) + gi'*A^s*W.c ) / ( sum(abs(gi'*W.G)) + gi'*W.c );
    end
    alpha_o = max(alpha_i);
else
    disp('Error in alpha_min.m: Set must be in either H-Rep or G-Rep')
    alpha_o = [];
end
    