function [x] = conszonohalfspaceIntersection_LP(z,temp_Y)

% conszonohalfspaceIntersection_LP computes the intersection $z \cap temp_Y$ by solving a linear program

% Inputs: z (constrained zonotope) (in CG-Rep) (struct variable with c,G,A,b parameters) 
% c-center of zonotope, G-generator matrix,A,b satisfying 
% z = {c +G\xi, ||\xi||_{\infty} < 1, A\xi = b (constraints)}

% temp_Y = \{y \in temp_Y, Py \leq q \} (Halfspace - A set of hyperplanes)

% Returns the intersection $x = z \cap temp_Y$ 

x = z; 

temp_Y = temp_Y.H; % Halfspace
P = temp_Y(:,1:end-1);
q = temp_Y(:,end);

for j = 1:size(temp_Y,1) % For every hyperplane
    xi1_ = sdpvar(size(x.G,2),1);
    objs = -P(j,:)*(x.c + x.G * xi1_); % Finds a point on $z$ orthogonal to the hyperplane
    cons = [];
    cons = [cons, -1 <= xi1_ <= 1, x.A*xi1_ == x.b];        
    opts = sdpsettings('solver','gurobi','verbose',1);
    O1 = optimizer(cons,objs,opts,[],[xi1_]);
    sol = O1{[]}; % Solves the LP
    d_zono = P(j,:)*(x.c + x.G *sol);
     if (d_zono >= q(j)) % Cons zono intersects with hyperplane
      
        d_max = temp_Y(j,end)-temp_Y(j,1:end-1)*x.c+sum(abs(temp_Y(j,1:end-1)*x.G));
        
            x.A = [x.A zeros(size(x.A,1),1); temp_Y(j,1:end-1)*x.G d_max/2]; % Adds hyperplane constraints
            x.b = [x.b; temp_Y(j,end)-temp_Y(j,1:end-1)*x.c-d_max/2];
            x.c = x.c;
            x.G = [x.G, zeros(size(x.G,1),1)];
    end
end