function [x,tcalc_all] = conszonohalfspaceIntersection_LP_comp(z,Y)

% conszonohalfspaceIntersection_LP_comp computes the intersection $z \cap Y$ 
% by solving a linear program and associated time.

% Inputs: z (constrained zonotope) (in CG-Rep) (struct variable with c,G,A,b parameters) 
% c-center of zonotope, G-generator matrix, A, b satisfying 
% z = {c +G\xi, ||\xi||_{\infty} < 1, A\xi = b (constraints)}

% Y = \{y \in Y, Py \leq q \} (Halfspace - A set of hyperplanes)

% Returns the intersection $x = z \cap Y$ and the computation time $tcalc_all$.

x = z;

temp_Y = Y.H; % Halfspace
P = temp_Y(:,1:end-1);
q = temp_Y(:,end);
tcalc_LP_sum = 0;
n_runs = 100; % Averaging computational time over n_runs.

for j = 1:size(Y.H,1)
        for m = 1:n_runs
            x = z;
            xi1_ = sdpvar(size(x.G,2),1);
			objs = -P(j,:)*(x.c + x.G * xi1_); 
			cons = [];
			cons = [cons, -1 <= xi1_ <= 1, x.A*xi1_ == x.b];        
			opts = sdpsettings('solver','gurobi','verbose',1);
			O1 = optimizer(cons,objs,opts,[],[xi1_]);
			tstart = tic;
            sol = O1{[]};
            d_zono = 0;
            d_zono = P(j,:)*(x.c + x.G *sol);                          
            if (d_zono >= q(j)) % Cons zono intersects with hyperplane

                d_max = Y.H(j,end)-Y.H(j,1:end-1)*x.c+sum(abs(Y.H(j,1:end-1)*x.G));

                    x.A = [x.A zeros(size(x.A,1),1); Y.H(j,1:end-1)*x.G d_max/2];
                    x.b = [x.b; Y.H(j,end)-Y.H(j,1:end-1)*x.c-d_max/2];
                    x.c = x.c;
                    x.G = [x.G, zeros(size(x.G,1),1)];
            end
            tcalc_iter(m) = toc(tstart);
        end                          
            tcalc_LP_sum = tcalc_LP_sum + mean(tcalc_iter);
end
tcalc_all = tcalc_LP_sum;
end
