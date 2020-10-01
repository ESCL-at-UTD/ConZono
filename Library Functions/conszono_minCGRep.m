function [x_r] = conszono_minCGRep(r,maxIter)
% conszono_minCGRep computes the minimal CG Rep(removes redundant generators and constraints) of a constrained zonotope 

% Inputs - r (Constrained zonotope) in CG-Rep satisfies r = {c + G\xi, ||\xi||_{\infty} \leq 1, A\xi = b}
% maxIter - Maximum number of iterations for refinement of intervals R, E

% Returns a constrained zonotope (x_r) in minimal CG-Rep as a struct variable. 

[temp, ind] = rrefcp([r.A,r.b]); % Computes the reduced row echelon form of [A b].
r.A = temp(:,1:size(r.A,2));
r.b = temp(:,size(r.A,2)+1);
r.G = r.G(:,ind(1:end-1));

x_r = r;
nc = size(x_r.A,1); % Determines the number of constraints
ng = size(x_r.A,2); % Determines the number of generators
[R,E] = Bounds(x_r,maxIter); % Computes the bounds R,E for the generators
redund = []; 
for j = 1:ng % Increments through the generators
    if (R(j,1) >= -1) && (R(j,2) <= 1) % Checks for a redundant generator
        redund = [redund;j]; % Concatenates all the redundant generators
    end
end

while length(redund) > 0 % Iterates until redundant generators are removed
 
    [R_ind,E_ind] = Bounds_ind(x_r,maxIter); % Computes R_ind,E_ind bounds of the generators for each of the constraints separately

    Ej1 = zeros(ng,nc);
    temp_Rmax_nc = zeros(nc,1);
    for j =1:nc % Computes the maximum of the absolute value among R_lb and R_ub.
       temp_Rmax_nc(j) = max(abs(R_ind(redund(1),j)),abs(R_ind(redund(1), nc + j)));
    end
    [~, ind_consrem] = min(temp_Rmax_nc); % Stores the constraint index that has min(max(R_lb, R_ub)) among all the constraints
        if x_r.A(ind_consrem,redund(1))  >= eps % Checks if the element is nonzero
            Ej1(redund(1),ind_consrem) = 1;
            Lambda_G = x_r.G*Ej1/x_r.A(ind_consrem,redund(1));
            Lambda_A = x_r.A*Ej1/x_r.A(ind_consrem,redund(1));
            x_r.c = x_r.c + Lambda_G*x_r.b; % Lines 39-42 performs order reduction
            x_r.G = x_r.G - Lambda_G*x_r.A; 
            x_r.A = x_r.A - Lambda_A*x_r.A;
            x_r.b = x_r.b - Lambda_A*x_r.b;
            x_r.A(ind_consrem,:) = []; % Lines 43-46 removes removes 1 generator and 1 constraint.
            x_r.A(:,redund(1)) = [];
            x_r.b(ind_consrem,:) = [];
            x_r.G(:,redund(1)) = [];
        end
    
    nc = size(x_r.A,1);
    ng = size(x_r.A,2);

    [R,E] = Bounds(x_r,maxIter); % Recomputes the R and E bounds.
     redund = [];
     for j = 1:ng
         if (R(j,1) >= -1) && (R(j,2) <= 1) % Checks for a redundant generator
             redund = [redund;j]; % Concatenates all the redundant generators
         end
     end              
end

end