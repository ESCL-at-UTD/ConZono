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
r = conZono;
r.c = [5;5];
r.G = eye(2);

% Input set
u = conZono;
u.c = 0;
u.G = 1;

% State Constraint Set
x_lb = [-4.05; -2.75];
x_ub = [6; 6];
X = Polyhedron('lb',x_lb,'ub',x_ub);


x = conZono;
x.c = (x_ub+x_lb)/2;
x.G = diag((x_ub-x_lb)/2);

%%
% Number of steps in reachability analysis
N = 6;

%%
figure;hold on
plot(x,'k',0.2)
plot(r)
xlim([-8 8])
ylim([-8 8])
drawnow

for k = 1:N
    r = A*r + B*u;
    r = halfspaceIntersection(r,X);
    r.plot
%     ri = innerApprox(x,r,'Infinity Norm');
    ri = innerApproxPointContain(x,r,'Infinity Norm',r.c);
%     tic
%     ri = innerApprox(x,r,'Hausdorff');
%     toc
    plot(ri,'b',0.2)
    plot(r.c(1),r.c(2),'.k')
    
    % Add point containment requirement
    
    xlim([-8 8])
    ylim([-8 8])
    drawnow
end
