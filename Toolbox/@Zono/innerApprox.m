function out = innerApprox(obj,n_r)
obj0 = copy(obj);

n_g = obj0.nG;
Norms = vecnorm(obj0.G);
[~,indx_sort] = sort(Norms,'descend');
obj0.G = obj0.G(:,indx_sort);

Z1_prim.G = obj0.G(:,1:n_r);
Z2_aux.G = obj0.G(:,n_r+1:end);

% Compute magnitude of dot product between generators in Z1 and Z2
alpha_abs = zeros(n_r,n_g-n_r);
alpha = zeros(n_r,n_g-n_r);
for i = 1:n_r
    for j = 1:n_g-n_r
        alpha_abs(i,j) = abs(Z1_prim.G(:,i)'*Z2_aux.G(:,j));
        alpha(i,j) = Z1_prim.G(:,i)'*Z2_aux.G(:,j);
    end
end
% Normalize the dot product with respect to largest dot product
alpha_norm = alpha*diag(1./max(alpha_abs,[],1));

T2 = zeros(n_r,n_g-n_r);
T2(alpha_norm==1) = 1;
T2(alpha_norm==-1) = -1;

T = [eye(n_r);T2'];

out = copy(obj)
out.c = obj0.c;
out.G = obj0.G*T;
end