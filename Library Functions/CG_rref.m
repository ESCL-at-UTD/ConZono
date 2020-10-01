function [z] = CG_rref(x)
% Computes rref (Reduced row echelon form with full pivoting) for a constrained zonotope.

% Inputs - x (Constrained zonotope) in CG-Rep satisfies x = {c + G\xi, ||\xi||_{\infty} \leq 1, A\xi = b}

% Returns z - a constrained zonotope with [A b] matrix in rref.

[temp_Ab_rref, ind] = rrefcp([x.A,x.b]); % Computes reduced row echelon form with full pivoting
z.A = temp_Ab_rref(:,1:size(x.A,2));
z.b = temp_Ab_rref(:,size(x.A,2)+1);
z.G = x.G(:,ind(1:end-1));
z.c = x.c;
end