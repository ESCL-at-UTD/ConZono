function [x, rep] = conszonohalfspaceIntersection_intar(z,Y,maxIter)
% conszonohalfspaceIntersection_intar computes the intersection $z \cap Y$ using Algorithm 1 from 
% "Constrained zonotopes:A new tool for set-based estimation and fault detection", Automatica, Scott et.al

% Inputs: z (constrained zonotope) (in CG-Rep) (struct variable with c,G,A,b parameters) 
% c-center of zonotope, G-generator matrix,A,b satisfying 
% z = {c +G\xi, ||\xi||_{\infty} < 1, A\xi = b (constraints)}
% Y = \{y \in Y, Py \leq q \} (Halfspace - A set of hyperplanes)
% maxIter - Maximum number of iterations for refinement of intervals R, E of x

% Returns the intersection $x = z \cap Y$ and rep (number of iterations 
% needed for refinement of intervals R, E of x)

x = z;

temp_Y = -Y.H; % Opposite extending halfspace

for m = 1:size(temp_Y,1)

% Add constraints of opposite extending hyperplane $H_{+}$

        d_max = temp_Y(m,end)-temp_Y(m,1:end-1)*x.c+sum(abs(temp_Y(m,1:end-1)*x.G));

        x.A = [x.A zeros(size(x.A,1),1); temp_Y(m,1:end-1)*x.G d_max/2];
        x.b = [x.b; temp_Y(m,end)-temp_Y(m,1:end-1)*x.c-d_max/2];
        x.c = x.c;
        x.G = [x.G, zeros(size(x.G,1),1)];
                             
nc = size(x.A,1);
ng = size(x.A,2);

e_lb = -ones(ng,1);
e_ub =  ones(ng,1);

E_old = [e_lb e_ub];

for rep = 1:maxIter
    r_lb = -inf*ones(ng,1);
    r_ub =  inf*ones(ng,1);

            for i = 1:nc
                for j = 1:ng
                    if abs(x.A(i,j)) >= eps % Machine precision                    
                    rj_temp = [0 0];
                    for k = 1:ng
                        if k ~= j
                            if x.A(i,k)/x.A(i,j) >= 0
                                rj_temp(1) = rj_temp(1) + x.A(i,k)/x.A(i,j)*e_lb(k);
                                rj_temp(2) = rj_temp(2) + x.A(i,k)/x.A(i,j)*e_ub(k);
                            else
                                rj_temp(1) = rj_temp(1) + x.A(i,k)/x.A(i,j)*e_ub(k);
                                rj_temp(2) = rj_temp(2) + x.A(i,k)/x.A(i,j)*e_lb(k);
                            end
                        end
                    end                        
                    rj_temp = x.b(i)/x.A(i,j) - [rj_temp(2) rj_temp(1)];
                    r_lb(j) = max(r_lb(j),rj_temp(1));
                    r_ub(j) = min(r_ub(j),rj_temp(2));
                    e_lb(j) = max(e_lb(j),r_lb(j));
                    e_ub(j) = min(e_ub(j),r_ub(j));
                    end
                end
            end
            R = [r_lb r_ub];
            E = [e_lb e_ub];
            sum_feas = 0;
            for j = 1:ng
                if (e_lb(j) <= e_ub(j)) % Checks if E is feasible
                sum_feas = sum_feas + 1;
                end
            end
            if sum_feas ~= ng
                break
            end
            if isequal(E, E_old) % Checks if the interval set E has been refined since the last iteration       
                break
            else
                E_old = E;
            end
end
            x.G(:,end) = []; % Remove the new generator
            x.A(nc,:) = []; % Remove the added constraint
            x.A(:,end) = []; % Remove the constraints on the new generator
            x.b(nc) = []; % Remove the added constraint                        
            
			if (sum_feas == ng) % If $z \cap H_{+} \neq \emptyset$                                         
                    % Add new generator and constraints corresponding to true
                    % hyperplane H_{-}
                        temp2_Y = -temp_Y(m,:); % Original halfspace
                        d_max = temp2_Y(1,end)-temp2_Y(1,1:end-1)*x.c+sum(abs(temp2_Y(1,1:end-1)*x.G));

                        x.A = [x.A zeros(size(x.A,1),1); temp2_Y(1,1:end-1)*x.G d_max/2];
                        x.b = [x.b; temp2_Y(1,end)-temp2_Y(1,1:end-1)*x.c-d_max/2];
                        x.c = x.c;
                        x.G = [x.G, zeros(size(x.G,1),1)];                
            end                
end                

end