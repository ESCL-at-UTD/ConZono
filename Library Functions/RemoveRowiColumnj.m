function [x_r] = RemoveRowiColumnj(x,i,j)
% RemoveRowiColumnj removes the ith row and the jth column of a constrained zonotope.

% Inputs - x (Constrained zonotope) in CG-Rep satisfies x = {c + G\xi, ||\xi||_{\infty} \leq 1, A\xi = b}
%          i - index of ith row to be removed
%          j - index of jth column to be removed
% Returns x_r - a constrained zonotope (x_r) with one less constraint and
% generator.

nc = size(x.A,1);
ng = size(x.A,2);

Eji = zeros(ng,nc);
Eji(j,i) = 1;
Lambda_G = x.G*Eji/x.A(i,j);
Lambda_A = x.A*Eji/x.A(i,j);
x_r.c = x.c + Lambda_G*x.b; %Lines 17-20 performs reduction
x_r.G = x.G - Lambda_G*x.A;
x_r.A = x.A - Lambda_A*x.A;
x_r.b = x.b - Lambda_A*x.b;
x_r.A(i,:) = []; % Lines 21-24 removes 1 generator and 1 constraint
x_r.A(:,j) = [];
x_r.b(i,:) = [];
x_r.G(:,j) = [];

end