%% Example System
Ac = [-1 -4; 4 -1];
Bc = [1;0];
sys_c = ss(Ac,Bc,eye(2),zeros(2,1));
dt = 0.1;
sys_d = c2d(sys_c,dt);
A = sys_d.A;
B = sys_d.B;

%% 
% Initial set for reachability
X0 = conZono;
X0.c = [5;5];
X0.G = eye(2);

% Input set
U = conZono;
U.c = 0;
U.G = 1;

% State Constraint Set
x_lb = [-4.05; -2.75];
x_ub = [6; 6];
X = Polyhedron('lb',x_lb,'ub',x_ub);

x = conZono;
x.c = (x_ub+x_lb)/2;
x.G = diag((x_ub-x_lb)/2);

%% Computes reachable sets
% Number of steps in reachability analysis
N = 6;

% compute reachable sets without intersection with feasible set
R = reachDT(A,B,X0,U,N);

% compute reachable sets with intersection with feasible set X
R_intersect = reachDT(A,B,X0,U,N,X);

%% Plots
figure;hold on;
plot(x,'k',0.2);
plot(X0);
xlim([-8 8]);
ylim([-8 8]);
drawnow

for k = 1:N
    R{k+1}.plot %reachable set
    R_intersect{k+1}.plot('b',1.0); %reachable set with intersection
    plot(R{k+1}.c(1),R{k+1}.c(2),'.k'); %center of the reachable set
    
    xlim([-8 8])
    ylim([-8 8])
    drawnow
end
