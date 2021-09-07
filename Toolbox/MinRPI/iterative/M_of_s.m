function [M_s] = M_of_s(W,A,s,MAX_ITER)
    if ~isprop(W,'H') && ~isfield(W,'G') && ~isprop(W,'G')
        error('Set must be in either H-Rep or G-Rep');
    end
    
    persistent n A_pow_i I sum_plus_j sum_minus_j
    
    %initialize persistent variables
    if s==1
        n = size(A,1);
        
        A_pow_i = cell(MAX_ITER,1);
        A_pow_i{1} = A;
        
        I = eye(n);
        
        sum_plus_j = zeros(n,1);
        sum_minus_j = zeros(n,1);
    else %create new value
        A_pow_i{s} = A^s;       
    end
    
    
    
    sum_max = zeros(n,1);
    for j = 1:n      
        a = A_pow_i{s}'*I(:,j);
        if isprop(W,'H')
            hW_a_plus = support(W,a);
            hW_a_minus = support(W,-a);
        else %G is a field of W
            hW_a_plus = sum(abs(a'*W.G)) + a'*W.c;
            hW_a_minus = sum(abs(a'*W.G)) - a'*W.c;
        end
        sum_plus_j(j) = sum_plus_j(j) + hW_a_plus;
        sum_minus_j(j) = sum_minus_j(j) + hW_a_minus;
        
        sum_max(j) = max(sum_plus_j(j),sum_minus_j(j));
    end
    M_s = max(sum_max);
    
end
