%% Application to Reachability Analysis 

close all;
clear all;

%% Plot Settings 
set(0,'defaultLineLineWidth', 2)
set(0,'defaultAxesFontName' , 'Times')
set(0,'defaultTextFontName' , 'Times')
set(0,'defaultAxesFontSize' , 14)
set(0,'defaultTextFontSize' , 14)
set(0,'defaulttextinterpreter','latex')
set(0,'defaultAxesGridLineStyle','-.')

%% System Definition

close all;
clear all;

% System size
Sys.Nx = 3;  % Number of states (Position, Velocity, Remaining Energy)
Sys.Nu = 3;  % Number of inputs (Acceleration, Braking, Load)
Sys.Ny = 6;  % Number of outputs (All States, All Inputs)
Sys.Nr = 4;  % Number of references (Position, Acceleration, Braking, Load)

% System parameters
Sys.a = 1;   % Acceleration coefficient
Sys.b = 1;   % Deceleration coefficient
Sys.M = 1;   % Mass of the system
Sys.dt = 1;  % Time step

% System dynamics
Sys.A = [1 Sys.dt 0; 0 1 0; 0 0 1];
Sys.B = [0 0 0; Sys.a*Sys.dt/Sys.M -Sys.b*Sys.dt/Sys.M 0; -1 -1 -1];

% System outputs y = [x1;x2;x3;u1;u2;u3];
Sys.C = [eye(Sys.Nx);zeros(Sys.Nu,Sys.Nx)];
Sys.D = [zeros(Sys.Nx,Sys.Nu);eye(Sys.Nu)];

% System reference outputs r = [x1;u1;u2;u3];
Sys.E = [1 0 0;zeros(Sys.Nu,Sys.Nx)];
Sys.F = [zeros(1,Sys.Nu);eye(Sys.Nu,Sys.Nu)];

% System Initial Condition
Sys.x0 = [0;0;100];

% Output constraints
Sys.Y = Polyhedron('lb',[-1 -20 0 0 0 0],'ub',[105 20 Sys.x0(3) 1 1 1]);
Sys.P = Sys.Y.H(:,1:end-1);
Sys.q = Sys.Y.H(:,end);

% State constraints 
Sys.X = projection(Sys.Y,1:Sys.Nx);

% Input constraints 
Sys.U = projection(Sys.Y,Sys.Nx+1:Sys.Nx+Sys.Nu);

% Terminal constraint
Sys.dx = 1;
Sys.dv = 1;
Sys.Xf = Polyhedron('lb',[-Sys.dx -Sys.dv 0],'ub',[Sys.dx Sys.dv Sys.x0(3)]);
Sys.Pf = Sys.Xf.H(:,1:end-1);
Sys.qf = Sys.Xf.H(:,end);
% 
%% Centralized MPC
Cent = Sys;
Cent.Name = 'Cent';
Cent.N = 100;  % horizon
Cent.weight = [1e2; 1e0; 1e0; 1e2];

Cent.rdes = zeros(Cent.Nr,Cent.N);
Cent.rdes(1,10:9*Cent.N/10) = 100;  % position reference
Cent.rdes(4,31:70) = 1;  % load reference

%% Upper Level System Structure

Up = Cent;
Up.dt = 10;
Up.nu = Up.dt/Cent.dt;

Up.A = Sys.A^Up.nu;
Up.B = Sys.B;
for i = 1:Up.nu-1
    Up.B = Up.B + Sys.A^i*Sys.B;
end

Sys.Ys = Sys.Y;

% Robustness to inter-sample behavior
for i = Up.nu:-1:1
    R = 0;
    for j = 0:i-1;
        R = R + Sys.A^j*Sys.B;
    end
    Pi = Sys.P*[Sys.C*Sys.A^i Sys.C*R + Sys.D];
    Yi = Polyhedron([Pi;Sys.P],[Sys.q;Sys.q]);
    Sys.Ys = intersect(Yi,Sys.Ys);
    Sys.Ys = minHRep(Sys.Ys);
end
%
Up.P = Sys.Ys.H(:,1:end-1);
Up.q = Sys.Ys.H(:,end);

Up.Name = 'Up';
Up.N = Cent.N/Up.nu;  % horizon
Up.weight = Cent.weight;

%% Upper level MPC

x_ = sdpvar(repmat(Up.Nx,1,Up.N+1),ones(1,Up.N+1));
u_ = sdpvar(repmat(Up.Nu,1,Up.N+1),ones(1,Up.N+1));
y_ = sdpvar(repmat(Up.Ny,1,Up.N),ones(1,Up.N));
r_ = sdpvar(repmat(Up.Nr,1,Up.N),ones(1,Up.N));
rdes_ = sdpvar(repmat(Up.Nr,1,Up.N),ones(1,Up.N));
t_ = sdpvar(repmat(1,1,Up.N+1),ones(1,Up.N+1));
tf_ = sdpvar(repmat(1,1,Up.N+1),ones(1,Up.N+1));

objs = 0;
cons = [];
for k = 1:Up.N
    for i = 1:Up.Nr
        objs = objs + Up.weight(i)*norm((rdes_{k}(i)-r_{k}(i))*(1-t_{k}),2)^2;
    end
    cons = [cons, x_{k+1} == Up.A*x_{k} + Up.B*u_{k}];
    cons = [cons, y_{k} == Up.C*x_{k} + Up.D*u_{k}];
    cons = [cons, r_{k} == Up.E*x_{k} + Up.F*u_{k}];
    cons = [cons, Up.P*y_{k} <= Up.q + ones(length(Up.q),1)*1e3*t_{k}];
    cons = [cons, Up.Pf*x_{k+1} <= Up.qf + ones(length(Up.qf),1)*1e3*tf_{k+1}];
end

opts = sdpsettings('solver','gurobi','gurobi.BarConvTol',1e-9,'gurobi.OptimalityTol',1e-9);  % Solve with gurobi
Up.Controller = optimizer(cons,objs,opts,{x_{1},rdes_{:},t_{:},tf_{:}},[x_,u_]); % Create exported optimization model object

rng(1)
Up.x0 = Cent.x0;
Up.rdes = mean(reshape(Cent.rdes,Up.Nr,Up.nu,Up.N),2);
Up.rdes = squeeze(Up.rdes);
Up.t = zeros(1,Up.N+1);
Up.t(end) = 1;
Up.tf = ones(1,Up.N+1);
Up.tf(end) = 0;

[ Up ]   = Call_Controller_Up( Up, 0, 0 );

%% Wayset computation (Backward reachability set)
% Method 1 - Zonotope Halfspace Intersection (Results in Table 2)

A_rev = inv(Sys.A);
B_rev = -inv(Sys.A)*Sys.B; 

u.ub = [1;1;1];
u.lb = [0;0;0];
x.ub = [105;20;105];
x.lb = [-1;-20;0];

u.c = (u.ub+u.lb)/2;
u.G = diag((u.ub-u.lb)/2);
u.A = zeros(0,size(u.G,2));
u.b = [];

x.c = (x.ub+x.lb)/2;
x.G = diag((x.ub-x.lb)/2);
x.A = zeros(0,size(x.G,2));
x.b = [];

k1 = 3;
N = Up.dt;

tcalc_ZH_wyset = []
n_runs = 1;
for j = 1:n_runs
    tstart = tic;
    r_ZH.c = Up.x(:,k1 + 3);
    r_ZH.G = zeros(Sys.Nx,0);
    r_ZH.A = [];
    r_ZH.b = [];

    for i = 1:N
        r_ZH.c = A_rev*r_ZH.c + B_rev*u.c;
        r_ZH.G = [A_rev*r_ZH.G B_rev*u.G];
        r_ZH.A = [r_ZH.A zeros(size(r_ZH.A,1),size(B_rev*u.G,2))];

        % State set intersection
        [r_ZH] = halfspaceIntersection(r_ZH,Sys.X);   
    end
    tcalc_ZH_wyset(j) = toc(tstart);
end

tcalc_ZH_wyset_mean = mean(tcalc_ZH_wyset)

size_ZH_A = size(r_ZH.A) % [7 37]

Wyset_HRep_orig = conszono_minHRep(r_ZH,1)

% Wyset_HRep_orig.volume

% Redundancy check 
n_runs = 1;
maxIter = 100;

tstart = [];
for j = 1:n_runs
    tstart = tic;
    % Redundancy Check 
    [x_r_ZH] = conszono_minCGRep(r_ZH,maxIter);
    tcalc_ZH_red(j) = toc(tstart);
end
tcalc_ZH_red_mean = mean(tcalc_ZH_red)

[zh_size] = size(r_ZH.A)
[zh_red_size] = size(x_r_ZH.A)

%% Method 2 - Generalized Intersection 

n_runs = 1;
for j = 1:1:n_runs
    tstart(j) = tic;
    r_GI.c = Up.x(:,k1 + 3);
    r_GI.G = zeros(Sys.Nx,0);
    r_GI.A = [];
    r_GI.b = [];
    
    r_GI_arr.c = [];
    r_GI_arr.G = {};
    r_GI_arr.A = {};
    r_GI_arr.b = {};
    
    r_GI_arr.c(:,Up.nu+1) = r_GI.c;
    r_GI_arr.G(Up.nu+1) = mat2cell(zeros(0,0),0,0);
    r_GI_arr.A(Up.nu+1) = mat2cell(zeros(0,0),0,0);
    r_GI_arr.b(Up.nu+1) = mat2cell([],0,0);
    
    for i = 1:N
        r_GI.c = A_rev*r_GI.c + B_rev*u.c;
        r_GI.G = [A_rev*r_GI.G B_rev*u.G];
        r_GI.A = [r_GI.A zeros(size(r_GI.A,1),size(B_rev*u.G,2))];

        % State set intersection
        [r_GI] = generalizedIntersection(r_GI,x,eye(Sys.Nx));
%         [r_GI] = conszono_minCGRep(r_GI,maxIter);
        r_GI_arr.c(:,i) = r_GI.c;
        r_GI_arr.G(i) = mat2cell(r_GI.G,size(r_GI.G,1),size(r_GI.G,2));
        r_GI_arr.A(i) = mat2cell(r_GI.A,size(r_GI.A,1),size(r_GI.A,2));
        r_GI_arr.b(i) = mat2cell(r_GI.b,size(r_GI.b,1),1);    
    end
    tcalc_GI_wyset(j) = toc(tstart(j));
end

tcalc_GI_wyset_mean = mean(tcalc_GI_wyset)

size_A = size(r_GI.A)

% r_GI_HRep = conszono_minHRep(r_GI,1);
% r_GI_HRep.volume;

% Redundancy check 

n_runs = 1;

for j = 1:n_runs
tstart = [];
tstart = tic;
[x_r_GI] = conszono_minCGRep(r_GI,maxIter);
tcalc_GI_red(j) = toc(tstart);
end

tcalc_GI_red_mean = mean(tcalc_GI_red)
[gen_size] = size(r_GI.A)
[gen_red_size] = size(x_r_GI.A)

%% Method 3 - Constrained Zonotope Halfspace Intersection using Linear Program

tcalc_LP_wyset = [];
n_runs = 1;
n_runs2 = 1;
for j = 1:1:n_runs
    r_LP = [];
    r_LP.c = Up.x(:,k1 + 3);
    r_LP.G = zeros(Sys.Nx,0);
    r_LP.A = [];
    r_LP.b = [];
%    tcalc_rchset_all = 0;
    for i = 1:N        
%        tstart_rchset = tic;
        for k = 1: n_runs2
        r_LP2.c = A_rev*r_LP.c + B_rev*u.c;
        r_LP2.G = [A_rev*r_LP.G B_rev*u.G];
        r_LP2.A = [r_LP.A zeros(size(r_LP.A,1),size(B_rev*u.G,2))];
        end
%        tcalc_rchset_all = tcalc_rchset_all + (toc(tstart_rchset)/n_runs); 
        r_LP.c = r_LP2.c;
        r_LP.G = r_LP2.G;
        r_LP.A = r_LP2.A;
        % State set intersection
        [r_LP] = conszonohalfspaceIntersection_LP(r_LP,Sys.X);   
%     [r_LP,tcalc_all(i)] = conszonohalfspaceIntersection_LP_comp(r_LP,Sys.X);
    end
end
% tcalc_all_sum = sum(tcalc_all) + tcalc_rchset_all

size_LP_A = size(r_LP.A) % [7 37]

% r_LP_HRep = conszono_minHRep(r_LP,1)
% r_LP_HRep.volume

%% Redundancy Detection 

n_runs = 1;
tstart = [];
% Redundancy Check 
    for k = 1:1:n_runs
        tstart = tic;
    [x_r_LP] = conszono_minCGRep(r_LP,maxIter);
    tcalc_lp_red(k) = toc(tstart);
    end
    
tcalc_lp_red_mean = mean(tcalc_lp_red)

[lp_size] = size(r_LP.A)
[lp_red_size] = size(x_r_LP.A)

%% Method 4 - Interval Arithmetic 

tcalc_IA_wyset = []
n_runs = 1;
maxIter = 100;
for j = 1:n_runs
    tstart = tic;
    r_IA.c = Up.x(:,k1+2+1);
    r_IA.G = zeros(3,0);
    r_IA.A = zeros(0,size(r_IA.G,2));
    r_IA.b = zeros(0,1);
    for i = 1:N
        r_IA.c = A_rev*r_IA.c + B_rev*u.c;
        r_IA.G = [A_rev*r_IA.G B_rev*u.G];
        r_IA.A = [r_IA.A zeros(size(r_IA.A,1),size(B_rev*u.G,2))];
        [r_IA,rep(i)] = conszonohalfspaceIntersection_intar(r_IA,Sys.X,maxIter);   
%         if size(r_IA.A,1) == 0 
%             [r_IA] = halfspaceIntersection(r_IA,Sys.X);       
%         else            
%         % State set intersection
%             [r_IA,rep(i)] = conszonohalfspaceIntersection_intar(r_IA,Sys.X,maxIter);   
%         end
    end
tcalc_IA_wyset(j) = toc(tstart);
end

rep
tcalc_IA_wyset_mean = mean(tcalc_IA_wyset)

size_IA_A = size(r_IA.A)
% 
% r_IA_HRep = conszono_minHRep(r_IA,1)
% r_IA_HRep.volume

%% Redundancy check 
n_runs = 1;
maxIter = 100;

tstart = [];
for j = 1:n_runs
    tstart = tic;
    % Redundancy Check 
    [x_r_IA] = conszono_minCGRep(r_IA,maxIter);
    tcalc_IA_red(j) = toc(tstart);
end
tcalc_IA_red_mean = mean(tcalc_IA_red)

[ia_size] = size(r_IA.A)
[ia_red_size] = size(x_r_IA.A)

%% Method 5 - H-Rep

t_calc_hrep = [];
size_hrep = [];
A = inv(Sys.A);
B = -inv(Sys.A)*Sys.B;
U = projection(Sys.Y,Sys.Nx+1:Sys.Nx+Sys.Nu);
X = projection(Sys.Y,1:Sys.Nx);
x3_hier = Up.x(:,3:11);

for k = 4
    tstart = tic;
    X_prime = x3_hier(:,k);

    for i = 1:Up.nu
        if i == 1
            X_prime = plus(A*X_prime,affineMap(U,B));
        else
            X_prime = plus(affineMap(X_prime,A),affineMap(U,B));
        end
        X_prime = intersect(X_prime,X);
%         X_prime = X_prime.minHRep;
        fprintf('%d ', i)
    end
    fprintf('\n')
    t_calc = toc(tstart);
    t_calc_hrep = [t_calc_hrep, t_calc];
    size_hrep = [size_hrep, size(X_prime.H,1)];
    
end
t_calc_hrep
size_hrep
mean(t_calc_hrep)

%% Redundancy detection (min H-Rep)
t_calc_minhrep = [];
size_minhrep = [];
A = inv(Sys.A);
B = -inv(Sys.A)*Sys.B;
U = projection(Sys.Y,Sys.Nx+1:Sys.Nx+Sys.Nu);
X = projection(Sys.Y,1:Sys.Nx);
tcalc_red = 0;

for k = 4
    tstart = tic;
    X_prime = x3_hier(:,k);

    for i = 1:Up.nu
        if i == 1
            X_prime = plus(A*X_prime,affineMap(U,B));
        else
            X_prime = plus(affineMap(X_prime,A),affineMap(U,B));
        end
        X_prime = intersect(X_prime,X);
%         X_prime = X_prime.minHRep;        
        fprintf('%d ', i)
    end
    fprintf('\n')
    tstart_red = tic;
    X_prime = X_prime.minHRep;        
    tcalc_red = toc(tstart_red);
    t_calc = toc(tstart);
    t_calc_minhrep = [t_calc_minhrep, t_calc];
    size_minhrep = [size_minhrep, size(X_prime.H,1)];
    
end
tcalc_red
t_calc_minhrep
size_minhrep
mean(t_calc_minhrep)

%% Evolution of wayset (Fig. 10)

load Aut_Zonotopes_Wayset_evolution_S13.mat 

figure('Position', [10 10 600 400]);

for i = Up.nu + 1:-1:1
      
   if i == Up.nu+1
%    r_S12.c = r_GI_arr.c(1:2,i);  
%    p_11 = plot3(i,r_S12.c(1,1),r_S12.c(2,1),'g.','MarkerSize',20);
   r_S12.c = r_GI_arr.c([1 3],i);  
   p_11 = plot3(i+4*Up.nu-1,r_S12.c(1,1),r_S12.c(2,1),'k.','MarkerSize',20);    
    hold on;
   elseif (i < Up.nu+1 && i > 1)
   r_S_all.G = cell2mat(r_GI_arr.G(i));
   r_S_all.A = cell2mat(r_GI_arr.A(i));
   r_S_all.b = cell2mat(r_GI_arr.b(i));

%    r_S12.c = [i;r_GI_arr.c(1:2,i)];
%    r_S12.G = [zeros(1,size(r_S_all.G,2));r_S_all.G(1:2,:)];
   r_S12.c = [i+4*Up.nu-1;r_GI_arr.c([1 3],i)];
   r_S12.G = [zeros(1,size(r_S_all.G,2));r_S_all.G([1 3],:)];
    r_S12.A = r_S_all.A;
   r_S12.b = r_S_all.b;
   r_evol(i) = conszono_minHRep(r_S12,1);
   p_10 = plot(r_evol(i),'color','g','alpha',0.2);    
   elseif i == 1
   r_S_all.G = cell2mat(r_GI_arr.G(i));
   r_S_all.A = cell2mat(r_GI_arr.A(i));
   r_S_all.b = cell2mat(r_GI_arr.b(i));   
% 
   r_S12.c = [i+4*Up.nu-1;r_GI_arr.c([1 3],i)];
   r_S12.G = [zeros(1,size(r_S_all.G,2));r_S_all.G([1 3],:)];
   r_S12.A = r_S_all.A;
   r_S12.b = r_S_all.b;
   r_evol(i) = conszono_minHRep(r_S12,1);
   p_1 = plot(r_evol(i),'color','r','alpha',0.2);           
   r_S12_opt.c = Up.x([1 3],5);  
   p_1_opt = plot3(i+4*Up.nu-1,r_S12_opt.c(1),r_S12_opt.c(2),'r.','MarkerSize',20);    
   end
   
end

% legend([p_11 p_10 p_1 p_1_opt],'$\mathbf{x}^{*}$','$Z_c(k + j), \forall j \in (1, N -1)$','$Z_c(k)$','$\mathbf{x}_{-}^{*}$','interpreter','Latex','NumColumns',4,'Location','northeast');
legend([p_11 p_10 p_1 p_1_opt],'$\mathbf{x}^{*}$','$Z_c(k + j), \forall j \in \{1, \cdots, N -1\}$','$Z_c(k)$','$\mathbf{x}_{-}^{*}$','interpreter','Latex','NumColumns',4,'Location','northoutside');

% [caz,cel] = view

caz = 38.8168;
cel = 15.0140;

view(caz,cel)

box on 
grid off
xlabel('Time [sec]');
ylabel('Position [m/sec]');
zlabel('Energy [kJ]');
set(gcf,'color','w');

% export_fig Wayset_S13_evol_time.pdf

%% Inner Approximation by a Box (Fig. 11 - Top Row)
Wyset_HRep_orig = conszono_minHRep(x_r_GI,2)

r_in_init.G = eye(Up.Nx);
r_in_init.c = x.c;
r_in_init.A = zeros(0,size(r_in_init.G,2));
r_in_init.b = [];

% Box 
b = r_in_init;

n_runs = 1;
for i = 1: n_runs
tstart_Box = tic;
[b_s] = ConZono_Containment_Opt_normval(b,x_r_GI,[],'inf');
tcalc_Box(i) = toc(tstart_Box);
end
tcalc_Box_mean = mean(tcalc_Box)

b_s.Box = Polyhedron('lb',-ones(size(b_s.G,2),1),'ub',ones(size(b_s.G,2),1));
Bs_set = plus(b_s.c,affineMap(b_s.Box,b_s.G));

R_vol_box_1step = (Bs_set.volume/Wyset_HRep_orig.volume)^(1/Sys.Nx) % 0.3517

Wyset_HRep_orig.contains(Bs_set)
Bs_set.contains(Up.x(:,k1+2))

%% Convex hull operation (CH(B \cup x_1^*(k_1 + 1))) (Fig. 11 Top Row)

x_Up.c = [Up.x(1,k1+1+1),Up.x(2,k1+1+1),Up.x(3,k1+1+1)]';
x_Up.G = zeros(Up.Nx,0);
x_Up.A = zeros(0,0);
x_Up.b = [];

for i = 1: n_runs
tstart_cvxhull = tic;
z_cvxhull_conszono_eqcons = cvxhull(b_s,x_Up);
tcalc_cvxhull(i) = toc(tstart_cvxhull);
end

tcalc_cvxhull_mean = mean(tcalc_cvxhull)

tcalc_box_nd_cvxhull_mean = tcalc_Box_mean + tcalc_cvxhull_mean;

Z_cvxhull_box_HRep = conszono_minHRep(z_cvxhull_conszono_eqcons,2);

str1 = ['Computation time to determine Box and convex hull operation are ',num2str(tcalc_Box_mean),' and ',num2str(tcalc_cvxhull_mean),' respectively ']; 
disp(str1)
str2 = ['Overall computation time is ',num2str(tcalc_box_nd_cvxhull_mean)]; 
disp(str2)

figure;
 plot(Wyset_HRep_orig,'color','g','alpha',0.2);
hold on;
 plot(Z_cvxhull_box_HRep,'color','k','alpha',0.5)
 hold on;
 plot3(Up.x(1,k1+1+1),Up.x(2,k1+1+1),Up.x(3,k1+1+1),'r.','MarkerSize',20);

size(z_cvxhull_conszono_eqcons.A) 
size(z_cvxhull_conszono_eqcons.b)

R_v_cvxhull_inapprox = (Z_cvxhull_box_HRep.volume/Wyset_HRep_orig.volume)^(1/Sys.Nx) %  

Z_cvxhull_box_HRep.contains(Up.x(:,k1+1+1))

%% Inner Approximation by a Zonotope containing the Up optimal solution (Fig. 11 Bottom Row)
% p = \infty

for i = 1:n_runs
    tstart_box_ptcontain = tic;
    [r_in_zon_norm_inf] = ConZono_Containment_Opt_normval(r_in_init,x_r_GI,Up.x(:,k1+2),'inf'); 
    tcalc_box_ptcontain(i) = toc(tstart_box_ptcontain); 
end
tcalc_box_ptcontain_mean = mean(tcalc_box_ptcontain)

disp(tcalc_box_ptcontain_mean)
r_in_zon_norm_inf_HRep = zono_minHRep(r_in_zon_norm_inf);

r_in_zon_norm_inf_HRep.contains(Up.x(:,k1+2))

R_v_zon_norm_inf_inapprox = (r_in_zon_norm_inf_HRep.volume/Wyset_HRep_orig.volume)^(1/Sys.Nx); % 0.3

Wyset_HRep_orig.contains(r_in_zon_norm_inf_HRep)

%% Plotting convexhull, Inner approximating zonotope (Up optimal solution enforcement) (Fig. 11 Bottom Row)
...Projection on to position-energy statespace

caz = [-39.1778];
cel = [7.1333];
pos_adj = 0.04;

h(1) = subplot(2,2,1);
h1_pos = get(h(1),'position')
h1_pos_upd = [h1_pos(1) h1_pos(2) (h1_pos(3) + pos_adj) h1_pos(4)];
set(h(1), 'position', h1_pos_upd);
p1 =  plot(Wyset_HRep_orig,'color','g','alpha',0.2);
hold on;
p2 = plot(Z_cvxhull_box_HRep,'color','k','alpha',0.5); 
p3 = plot(Bs_set,'color','m','alpha',0.8);drawnow;  
p4 = plot3(Up.x(1,k1+1+1),Up.x(2,k1+1+1),Up.x(3,k1+1+1),'r.','MarkerSize',20);
 grid off
 box on
set(gcf,'color','w');
 xlabel('Position [m]');   
 ylabel('Velocity [m/sec]');   
 zlabel('Energy [kJ]');   
 
view(caz,cel)

% Construct a legend with the data from the sub-plots
hL = legend([p1,p2,p3,p4],{'$Z_c$','$CH(B \cup \mathbf{x}_{-}^*)$','$B$','$\mathbf{x}_{-}^*$'},'interpreter','Latex');
% Positioning the legend
newPosition = [0.17 0.94 0.3 0.02];
newUnits = 'normalized';
hL.NumColumns =4;
hl.FontSize = 14;
hl.FontName = 'Times';

set(hL,'Position', newPosition,'Units', newUnits);

h(3) = subplot(2,2,3);
h3_pos = get(h(3),'position')
h3_pos_upd = [h3_pos(1) h3_pos(2) (h3_pos(3) + pos_adj) h3_pos(4)];
set(h(3), 'position', h3_pos_upd);
p5 = plot(Wyset_HRep_orig,'color','g','alpha',0.2); 
hold on;
p6 = plot(r_in_zon_norm_inf_HRep,'color','m','alpha',0.8);
 hold on;
 p7 = plot3(Up.x(1,k1+1+1),Up.x(2,k1+1+1),Up.x(3,k1+1+1),'r.','MarkerSize',20);
 grid off   
 box on
 set(gcf,'color','w');
 xlabel('Position [m]');   
 ylabel('Velocity [m/sec]');   
 zlabel('Energy [kJ]');   
view(caz,cel)

h(2) = subplot(2,2,2);
h2_pos = get(h(2),'position')
h2_pos_upd = [h2_pos(1) h2_pos(2) (h2_pos(3)-pos_adj) h2_pos(4)];
set(h(2), 'position', h2_pos_upd);

Wyset_HRep_orig_S12 = projection(Wyset_HRep_orig,[1:2]);
Wyset_HRep_orig_S13 = projection(Wyset_HRep_orig,[1 3]);
Z_cvxhull_HRep_S12 = projection(Z_cvxhull_box_HRep,1:2);
Z_cvxhull_HRep_S13 = projection(Z_cvxhull_box_HRep,[1 3]);
x_red_inapprox_HRep_S12 = projection(Bs_set,1:2);
x_red_inapprox_HRep_S13 = projection(Bs_set,[1 3]);

r_in_zon_norm_inf_HRep_S12 = projection(r_in_zon_norm_inf_HRep,[1 2]);

p(16) =  plot(Wyset_HRep_orig_S12,'color','g','alpha',0.2);
hold on;
p(17) = plot(Z_cvxhull_HRep_S12,'color','k','alpha',0.5); 
p(18) = plot(x_red_inapprox_HRep_S12,'color','m','alpha',0.8);drawnow
p(19) = plot(Up.x(1,k1+1+1),Up.x(2,k1+1+1),'r.','MarkerSize',20);
 grid off
 box on
 xlabel('Position [m]');   
 ylabel('Velocity [m/sec]');   
 set(gcf,'color','w');

h(4) = subplot(2,2,4);
h4_pos = get(h(4),'position')
h4_pos_upd = [h4_pos(1) h4_pos(2) (h4_pos(3) - pos_adj) h4_pos(4)];
set(h(4), 'position', h4_pos_upd);

p(16) =  plot(Wyset_HRep_orig_S12,'color','g','alpha',0.2); % 1
hold on;
p(17) = plot(r_in_zon_norm_inf_HRep_S12,'color','m','alpha',0.8);
p(19) = plot(Up.x(1,k1+1+1),Up.x(2,k1+1+1),'r.','MarkerSize',20);

grid off
 box on
 xlabel('Position [m]');   
 ylabel('Velocity [m/sec]');   
 set(gcf,'color','w');
 