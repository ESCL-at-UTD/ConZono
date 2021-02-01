function [out] = removerowicolumnj(obj,i,j)
out = copy(obj);

Eji = zeros(obj.nG,obj.nC);
Eji(j,i) = 1;
Lambda_G = obj.G*Eji/obj.A(i,j);
Lambda_A = obj.A*Eji/obj.A(i,j);
out.c = obj.c + Lambda_G*obj.b;
out.G = obj.G - Lambda_G*obj.A;
out.A = obj.A - Lambda_A*obj.A;
out.b = obj.b - Lambda_A*obj.b;
out.A(i,:) = [];
out.A(:,j) = [];
out.b(i,:) = [];
out.G(:,j) = [];

end