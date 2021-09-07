%% Original equations taken from  V. Raghuraman and J. P. Koeln, “Set operations and order reductions for constrained zonotopes.”, 2020.

function [RPIset] = RPI_1step(A,W,n_steps)

% computing candidate RPI set   
R = W;
for i = 1:n_steps
    R = R + A^i*W;
end
    
%% Create some of the main variables
n = size(A,1);
ng = R.nG;
nw = W.nG;
G = R.G;
Gw = W.G;
cw = W.c;

%% Identify sizes of each set of decision variables
n_gamma1 = ng*ng; %number of vars in Gamma1
n_gamma2 = ng*nw; %number of vars in Gamma2
n_phi = ng; % number of vars in phi
n_beta = ng; % number of vars in Beta
n_c = n; % number of var in c (center of the new set)
%% Constructing Equality constraints for the optimization problem

% Equation 28a)
%Block with columns of A*G:
AG = A*G;
H1_phi = sparse(n*ng,ng);
for i=1:ng
    H1_phi((i-1)*n+1:i*n,i) =  AG(:,i);
end


%Block with G's in the diagonal, related to Gamma1 variables + zeros for Gamma2, Beta and c variables:
H1 = [H1_phi, kron(speye(ng),-G), sparse(n*ng,n_gamma2+n_beta+n_c)];
b1 = sparse(size(H1,1),1);

%Equation 28b)
%Block with G's in the diagonal, related to Gamma2, and columns of Gw as
%the 'b' matrix
H2 = [sparse(n*nw,n_phi+n_gamma1),kron(speye(nw),G),sparse(n*nw,n_beta+n_c)];
b2 = reshape(Gw,[],1);

% Equation 28c)
H3 = [sparse(n,n_phi+n_gamma1+n_gamma2), -G, (eye(n)-A)];
b3 = cw;

%% Equation 28d)

%we need a new variable for each absolute value variable
n_gamma1_abs = n_gamma1;
n_gamma2_abs = n_gamma2;
n_beta_abs = n_beta;

%Here we construct the pattern for the sum of absolute variables for the
%first row of equation 28d) .
idx_gamma1 = zeros(1,n_gamma1_abs);
idx_gamma1(1:ng:end)=1;
idx_gamma2 = zeros(1,n_gamma2_abs);
idx_gamma2(1:ng:end)=1;
idx_beta = zeros(1,n_beta_abs);
idx_beta(1)=1;

H_abs = zeros(ng,n_gamma1_abs+n_gamma2_abs+n_beta_abs);
r_abs_part = [idx_gamma1 idx_gamma2 idx_beta];

%Now we assign it to the first row and then replicate that for the remaining rows
%by shifting the indexes 1position to the right each next row
for i=1:ng
    H_abs(i,:) =  r_abs_part;
    r_abs_part = circshift(r_abs_part,1);
end

%including the phi variables (which have no absolute value) from inequality
%28d), and zeros in place of the orginal variables
H_abs = [-eye(n_phi) sparse(ng,n_gamma1+n_gamma2+n_beta+n_c) sparse(H_abs)];

%Now we need to relate the original abs variables to their placeholders
%(t_i)
%thorugh inequalities (var_i - t_i <= 0) and (-var_i - t_i <= 0). That is,
%t_i = abs(var_i) and this is how we express that in a linear program.
n_t_vars =n_gamma1+n_gamma2+n_beta;
H_abs2 = sparse([sparse(n_t_vars,n_phi) +speye(n_t_vars) sparse(n_t_vars,n_c) -speye(n_t_vars);...
                 sparse(n_t_vars,n_phi) -speye(n_t_vars) sparse(n_t_vars,n_c) -speye(n_t_vars)]);

% PUtting these inequalities together
H = [H_abs ; H_abs2];
b = zeros(size(H,1),1);

% H = convert_abs(sparse(H_abs),H_phi); ==> this method also works. It does
% not require new t_i variables but requires addition of something of the
% order 2^(n_gamma1+n_gamma2+n_beta) inequalities, which quickly becomes
% prohibitive.
%% Putting equality constraints together
Heq = sparse([H1;H2;H3]);
Heq = [Heq, sparse(size(Heq,1),n_gamma1_abs+n_gamma2_abs+n_beta_abs)];
beq = [b1;b2;b3];

%% Include variable 't' (not the same t as before) to represent infinity norm on objective function
% We want to minimize the infinity norm of the phi vector. This is
% equivalent to minimizing t and adding constraints 1) -phi_i<=t and 2)
% pgi_i<=t.

%equality constraints are the same, just add a 0 for the new t variable
Heq = [Heq zeros(size(Heq,1),1)];

%we need two inequality rows for each phi_i:
Hobj = [ [speye(n_phi);-speye(n_phi)] sparse(2*n_phi,n_gamma1+n_gamma2+n_beta+n_c) sparse(2*n_phi,n_gamma1_abs+n_gamma2_abs+n_beta_abs) -ones(2*n_phi,1)];

%updates inequality constraints matrices
H = sparse([[H sparse(size(H,1),1)];Hobj]);
b = [b;zeros(2*n_phi,1)];

%% Objective function
%Final decision variables vector is given by x = [phi_1 ... phi_ng;
%Gamma1_c1_1.... Gamma1_c1_ng... Gamma1_cng_1... Gamma1_cng_ng;
% Gamma2_c1_1.... Gamma2_c1_ng... Gamma2_cnw_1... Gamma2_cnw_ng;
% Beta_1 ... Beta_ng ..... c_1....c_n;
% Gamma1_abs_c1_1.... Gamma1_abs_c1_ng... Gamma1_abs_cng_1... Gamma1_abs_cng_ng;
% Gamma2_abs_c1_1.... Gamma2_abs_c1_ng... Gamma2_abs_cnw_1... Gamma2_abs_cnw_ng;
% Beta_abs_1 ... Beta_abs_ng; t]
% Here the subscripts ci refer to the i-th column of a matrix.

%objective function
f = zeros(1,n_phi+n_gamma1+n_gamma2+n_beta+n_c+n_gamma1_abs+n_gamma2_abs+n_beta_abs+ 1);
f(end)=1; %minimize 't'


options = optimoptions(@linprog,'Display','off'); %linprog is just a dummy function here
[x,~,exitflag] = linprog_Gurobi(f,H,b,Heq,beq,[],[],options);
if exitflag == 1
    idx_1st_c = n_phi+n_gamma1+n_gamma2+n_beta+1;
    c = x(idx_1st_c:idx_1st_c+n_c-1);
    phi = x(1:n_phi);
else
    error("Error finding solution to the LP.");
end

G = G*diag(phi);

RPIset = conZono(c,G);
end

