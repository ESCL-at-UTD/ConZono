function out = pontryagin_diff_1step(obj1,obj2)

out = copy(obj1);

if obj1.n ~= obj2.n
    disp(['Cannot compute the Pontryagin difference between a zonotope in ',num2str(obj1.n),...
        ' dimensions and a zonotope in ',num2str(obj2.n),' dimensions!'])
else
    method=1;
    if method==1
        cd = sdpvar(obj1.n,1);
        S = sdpvar(obj1.nG+obj2.nG,1);
        Gamma = sdpvar(obj1.nG,obj1.nG+2*obj2.nG,'full');
        beta = sdpvar(obj1.nG,1);
        cons = [];
        cons = [cons, [[obj1.G obj2.G]*diag(S) obj2.G] == obj1.G*Gamma];
        cons = [cons, obj1.c - (cd+obj2.c) == obj1.G*beta];
        cons = [cons, abs(Gamma)*ones(obj1.nG+2*obj2.nG,1) + abs(beta) <= ones(obj1.nG,1)];
        cons = [cons, S >= 0];
        
        obj = norm([obj1.G obj2.G]*(1-S),'inf'); % Norm p = 1, 2, \infty
        opts = sdpsettings('solver','gurobi','verbose',0);
        
        Opt_solve = optimizer(cons,obj,opts,[],{cd,S}); % Create exported optimization model object
        [OUT_Opt,diagnostics] = Opt_solve([]); % Solves the LP
        out.c = cell2mat(OUT_Opt(1));
        S_val = cell2mat(OUT_Opt(2));
        
        Gout = [obj1.G obj2.G]*diag(S_val);
        %removing zero generators
        nG = size(Gout,2);
        cols_to_remove = zeros(nG,1);
        for k=1:nG
            if norm(Gout(:,k),Inf) ==0
                cols_to_remove(k) = k;
            end
        end
        Gout(:,cols_to_remove~=0)=[];
        
        out.G = Gout;
        
    else
        tic
        if false
            %final decision variables vector is x=[S;Gamma_c;Beta;cd;t]
            
            n = obj1.n;
            numel_Gamma = obj1.nG*(obj1.nG + 2*obj2.nG);
            numel_Beta = obj1.nG;
            num_cols_Gamma = (obj1.nG + 2*obj2.nG);
            numel_S = obj1.nG + obj2.nG;
            numel_cd = n;
            
            
            G1G2 = [obj1.G obj2.G];
            G1 = obj1.G;
            G2 = obj2.G;
            
            nc = (obj1.nG + obj2.nG);
            KG1G2  = zeros(obj1.n*nc,nc);
            
            for i=1:nc
                idx = (1 + (i-1)*n):(1 + (i-1)*n + n - 1);
                KG1G2(idx,i) = G1G2(:,i);
            end
            
            KGamma = kron(eye(nc),G1);
            
            Aeq1 = [KG1G2 , -KGamma];
            beq1 = [zeros(obj1.n*nc,1); reshape(-G2,[],1)];
            
            %equation of beta and cd
            %c1 -(cd + c2) = G1Beta ==> [G1 I][Beta;cd] = c1-c2
            Aeq2 = [G1 eye(n)];
            beq2 = obj1.c-obj2.c;
            
            Aeq = blkdiag(Aeq1,Aeq2,zeros(numel_Beta),0);
            beq = [beq1;beq2;zeros(numel_Beta,1),0];
            
            
            %Inequality |Gamma|1 + |Beta| <= 1
            
            %Inequality matrices after rearanging the inequality above to have
            %one variables column vector = [Beta;Gamma_c], where Gamma_c is the
            %matrix Gamma with all column vectors stacked above one another
            
            
            %Considering each row of the original matrix inequality, we have to
            %create 2 rows for each variable to get rid of the absolute values.
            %This results in all possible -1, +1 permutations of each variable
            % for each row
            P_row = double(permutations(1 + num_cols_Gamma));
            P_row(P_row==0) = -1;
            n_rows_P = size(P_row,1);
            
            A = zeros(obj1.nG*n_rows_P,numel_Beta + numel_Gamma);
            b = ones(obj1.nG*n_rows_P,1);
            
            %Now we have to "distribute" the permutation matrix above in the
            %corresponding rows and columns of the final A matrix. That is, the
            %permutation matrix is correct when each row of the original matrix
            %inequality is considered isolated. Now we have to relocate each
            %element of that in the corresponding rows and columns of A.
            gamma_col_indexes = numel_Beta+1:num_cols_Gamma:num_cols_Gamma;
            for i=1:obj1.nG
                
                col_indexes = [i (gamma_col_indexes + (i-1))]; %i corresponds to the Beta_i column
                
                for j=1:n_rows_P
                    %the first column corresponds to the Beta_i variable
                    A(((i-1)*n_rows_P) + j,col_indexes) = P_row(j,:);
                end
            end
            
            %lets include the other variables zeros in the matrix inequality
            Beta_block = A(:,1:numel_Beta);
            Gamma_block = A(:,numel_Beta+1:end);
            %here we have to invert the Beta and Gamma blocks since in the
            %final decision variable Gamma commes first
            Afinal = blkdiag(zeros(numel_S),[Gamma_block Beta_block], zeros(numel_cd));
            bfinal = [zeros(numel_S,1);b;zeros(numel_cd,1)];
            
            %Objective function: minimize(norm([G1 G2][1-S],inf)) ==>
            %[-G1G2 -1;G1G2 -1][S;t]<=[G1G2;-G1G2][1;1]
            Aobj = [-G1G2 zeros(n,numel_Gamma+numel_Beta+numel_cd) -ones(n,1); G1G2 zeros(n,numel_Gamma+numel_Beta+numel_cd) -ones(n,1)];
            bobj = [G1G2 ;-G1G2]*ones(2*n,1);
            
            Afinal = [Afinal;Aobj];
            bfinal = [bfinal;bobj];
            
            %we want to minimize t
            f = zeros(numel_S+numel_Gamma+numel_Beta+numel_cd+1,1);
            f(end) = 1;
            
            %final decision variables vector is x=[S;Gamma_c;Beta;cd;t]
            BigM = 1e10;
            lb = [zeros(numel_S,1);-ones(numel_Gamma,1);-ones(numel_Beta,1);-BigM*ones(numel_cd,1);-BigM];
            ub = [BigM*ones(numel_S,1);+ones(numel_Gamma,1);+ones(numel_Beta,1);+BigM*ones(numel_cd,1);+BigM];
        end
        
        
        cd = sdpvar(obj1.n,1);
        S = sdpvar(obj1.nG+obj2.nG,1);
        Gamma = sdpvar(obj1.nG,obj1.nG+2*obj2.nG,'full');
        beta = sdpvar(obj1.nG,1);
        cons = [];
        cons = [cons, [[obj1.G obj2.G]*diag(S) obj2.G] == obj1.G*Gamma];
        cons = [cons, obj1.c - (cd+obj2.c) == obj1.G*beta];
        cons = [cons, abs(Gamma)*ones(obj1.nG+2*obj2.nG,1) + abs(beta) <= ones(obj1.nG,1)];
        cons = [cons, S >= 0];
        
        obj = norm([obj1.G obj2.G]*(1-S),'inf'); % Norm p = 1, 2, \infty
        %         opts = sdpsettings('solver','gurobi','verbose',0);
        
        [Model,~,~,~] = export(cons, obj,sdpsettings('solver','linprog'));
        
        
        lp = Opt('f',Model.c,'Ae',Model.Aeq,'be',Model.beq,'A',Model.A,'b',Model.b,'lb',Model.lb,'ub',Model.ub); % Formulates the linear program
        
        sol = mpt_solve(lp); % Solves the LP.
        toc
        if sol.exitflag == 1
            cd = sol.xopt(1:obj1.n);
            S_val = sol.xopt(obj1.n+1:obj1.n+(obj1.nG +obj2.nG));
            
            out.c = cd;
            
            
            Gout = [obj1.G obj2.G]*diag(S_val);
            %removing zero generators
            nG = size(Gout,2);
            cols_to_remove = zeros(nG,1);
            for k=1:nG
                if norm(Gout(:,k),Inf) ==0
                    cols_to_remove(k) = k;
                end
            end
            Gout(:,cols_to_remove~=0)=[];
            
            out.G = Gout;
            
        else
            error("Error finding solution to the LP. i=%d",i);
        end
    end
    
    
    %     cd = sdpvar(obj1.n,1);
    % 	phi_ = sdpvar(obj1.nG+obj2.nG,1);
    % 	obj0 = copy(obj1);
    % 	obj0.c = cd + obj2.c;
    % 	obj0.G = [obj1.G obj2.G]*diag(phi_);
    %
    % 	cons = [phi_ >= 0];
    %
    % 	opts = sdpsettings('solver','gurobi','verbose',0);
    % 	objs = norm([obj1.G obj2.G]*(1-phi_),'inf');
    % 	[sol] = optimize(cons,objs,opts);
    % 	out.c = value(cd) + obj2.c;
    % 	out.G = [obj1.G obj2.G]*diag(value(phi_));
end
end


function P = permutations(N)
[P{1:N}] = ndgrid(logical([0 1]));
P = cat(N+1,P{N:-1:1});
P = reshape(P,[],N);
end
