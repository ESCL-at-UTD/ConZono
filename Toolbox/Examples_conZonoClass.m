%% Plot settings 
set(0,'defaultLineLineWidth', 2)
set(0,'defaultAxesFontName' , 'Times')
set(0,'defaultTextFontName' , 'Times')
set(0,'defaultAxesFontSize' , 18)
set(0,'defaultTextFontSize' , 18)
set(0,'defaulttextinterpreter','latex')
set(0,'defaultlegendinterpreter','latex')
set(0,'defaultAxesGridLineStyle','-.')

%% 1 - Constructing and plotting a constrained zonotope
z = conZono; % Creating a conZono object

z.c = [0;0];
z.G = [1 0 0; 0 1 0];
z.A = [1 1 1];
z.b = [1];

% Plotting a constrained zonotope
z.plot('b',0.1)
xlabel('$z_1$','interpreter','Latex');
ylabel('$z_2$','interpreter','Latex');
set(gcf,'color','w');
grid on;
box on;
legend('$Z_c$','interpreter','Latex');

%Alternative constructors
clearvars z;
c = [0;0];
G = [1 0 0; 0 1 0];
A = [1 1 1];
b = [1];

% c,g constructor
z = conZono(c,G);
figure; z.plot('b',0.1); hold on;
z.A = [1 1 1];
z.b = [1];
z.plot('r',0.1);

xlabel('$z_1$','interpreter','Latex');
ylabel('$z_2$','interpreter','Latex');
set(gcf,'color','w');
grid on;
box on;
legend('$Z_c$','interpreter','Latex');

% c,G,A,b constructor
clearvars z;
z = conZono(c,G,A,b);
figure; z.plot('b',0.1); hold on;
xlabel('$z_1$','interpreter','Latex');
ylabel('$z_2$','interpreter','Latex');
set(gcf,'color','w');
grid on;
box on;
legend('$Z_c$','interpreter','Latex');

%% 2 - Constructing and plotting the scalar multiplication of a zonotope
figure;hold on
z.plot
s = mtimes(z,2) % Scaling the constrained zonotope
% s = 2*z;
s.plot('b',0.1)
xlabel('$z_1$','interpreter','Latex');
ylabel('$z_2$','interpreter','Latex');
set(gcf,'color','w');
grid on;
box on;
legend('$Z_c$','$2*Z_c$','interpreter','Latex');

%% 3- Computing and plotting the Minkowski sum of two constrained zonotopes
x = copy(z);
y = copy(z);

figure;
x.plot('b',0.5)
hold on
w = x + y;
w.plot('m',0.1)
xlabel('$z_1$','interpreter','Latex');
ylabel('$z_2$','interpreter','Latex');
set(gcf,'color','w');
grid on;
box on;
legend('$\mathcal{X}, \mathcal{Y}$','$\mathcal{X} \oplus \mathcal{Y}$','interpreter','Latex');

%% 4- Computing and plotting Zonotope-hyperplane and constrained zonotope-hyperplane intersection 

x = conZono;
x.c = [0;0];
x.G = [1 1; 0 2];
x.A = zeros(0,size(x.G,2));
x.b = [];

% Generating Halfspace and Hyperplane
E = [3 1];
f = [3];
H_ = Polyhedron('H',[E f]);
H = Polyhedron('He',[E f]);

% Evaluating intersection
abs(f-E*x.c)
sum(abs(E*x.G))

% Computing Intersection of Zonotope and Halfspace
 
[z] = halfspaceIntersection(x,H_);

% Generating Second Halfspace and Hyperplane
E2 = [4 1];
f2 = [6];
H2_ = Polyhedron('H',[E2 f2]);
H2 = Polyhedron('He',[E2 f2]);

%  Plotting

H_mod = Polyhedron('H',[E f],'lb',[-3 -3],'ub',[3 3]);
H2_mod = Polyhedron('H',[E2 f2],'lb',[-3 -3],'ub',[3 3]);

figure('Position',[100 100 900 600])
subplot(1,2,1); hold on
p1 = plot(H_mod,'color',[0.9 0.7 1],'LineStyle','none');
x.plot('k',0.3);
p3 = plot([2 0],[-3 3],'color',[0.6 0.3 1],'LineWidth',3);
z.plot('c',0.7);
p5 = plot(x.c(1),x.c(2),'.k','markerSize',20);
set(gca, 'Layer','top')
plot([-2 -2.4],[-2 -0.8],'--','color',[0.6 0.3 1],'LineWidth',1);

annotation('doublearrow',[0.305 0.34],[0.523 0.545],'HeadStyle','plain','Head1Length',6,'Head1Width',6,'Head2Length',6,'Head2Width',6)
annotation('doublearrow',[0.289 0.17],[0.513 0.443],'HeadStyle','plain','Head1Length',6,'Head1Width',6,'Head2Length',6,'Head2Width',6)
annotation('textbox',[0.28 0.53 0.04 0.05],'string','$\bf{c}$','Interpreter','latex','lineStyle','none')
annotation('textbox',[0.31 0.488 0.04 0.05],'string','$d_1$','Interpreter','latex','lineStyle','none')
annotation('textbox',[0.227 0.44 0.04 0.05],'string','$d_2$','Interpreter','latex','lineStyle','none')

xlim([-3 3])
ylim([-3 3])
xlabel('$z_1$')
ylabel('$z_2$')
yticks(linspace(-2,2,3))

grid off
box on
axis square
leg = legend('$H_-$','$Z$','$H$','$Z_h$');
set(leg,'Interpreter','latex','location','northwest');

subplot(1,2,2); hold on
plot(H2_mod,'color',[0.9 0.7 1],'LineStyle','none');
x.plot('k',0.3);
plot([2.25 0.75],[-3 3],'color',[0.6 0.3 1],'LineWidth',3);
z.plot('c',0.7);
p5 = plot(x.c(1),x.c(2),'.k','markerSize',20);
set(gca, 'Layer','top')
annotation('textbox',[0.72 0.53 0.04 0.05],'string','$\bf{c}$','Interpreter','latex','lineStyle','none')

xlim([-3 3])
ylim([-3 3])
xlabel('$z_1$')
ylabel('$z_2$')
yticks(linspace(-2,2,3))

grid off
box on
axis square
leg = legend('$H_-$','$Z$','$H$','$Z_c$');

% leg = legend([p2 p4 p1 p3],'$Z$','$Z_c$','$H_-$','$H$');
set(leg,'Interpreter','latex','location','northwest');


%% 5- Plotting the generalized intersection of two zonotopes
y = copy(z);
y.b = 2;
y.c = [0.5;0];
figure;hold on
z.plot
y.plot('b',0.7)
x = generalizedIntersection(z,y,eye(2));
x.plot('k',0.7)

xlabel('$z_1$','interpreter','Latex');
ylabel('$z_2$','interpreter','Latex');
set(gcf,'color','w');
grid on;
box on;
legend('$\mathcal{Z}$','$\mathcal{Y}$','$\mathcal{Z} \cap \mathcal{Y}$','interpreter','Latex','Location','southwest');

%% 6- Convex hull 

% convex hull of zonotopes
x = conZono;

x.c = [0;0];
x.G = [0 1 0; 1 1 2];
x.A = zeros(0,3);
x.b = [];

y = conZono;
y.c = [-5;0];
y.G = 0.5*[-1 2 -4;1 1 3];
y.A = zeros(0,3);
y.b = [];

hull_uncons = convexHull(x,y);
size(hull_uncons.A)

% Check 

containCheck(x,hull_uncons)
containCheck(y,hull_uncons)

% Computing the convex hull of Constrained zonotopes
% Intersecting x with a halfspace 
E1 = [1 1];
f1 = [0];
H1 = Polyhedron('H',[E1 f1]);

[x_cons] = halfspaceIntersection(x,H1);

% Intersecting y with a halfspace 
E2 = [-2.5 1];
f2 = 9.5;
H2 = Polyhedron('H',[E2 f2]);

[y_cons] = halfspaceIntersection(y,H2);

% Computing the cvxhull
hull_cons = convexHull(x_cons,y_cons);

% Constraining the length of the halfspaces H1, H2, H3,H4 for easy plotting
H1_cons = conZono;
H1_cons.G = 1.25*[-1.5 1.5]' ;
H1_cons.c = [0;0];

Box_H1_cons = Polyhedron('lb',-ones(size(H1_cons.G,2),1),'ub',ones(size(H1_cons.G,2),1));
X_H1_cons = plus(H1_cons.c,affineMap(Box_H1_cons,H1_cons.G));

H2_cons = conZono;
H2_cons.c = [-4 -0.5]' ;
H2_cons.G = 1.25*[0.75 ; 1.25*1.5];

Box_H2_cons = Polyhedron('lb',-ones(size(H2_cons.G,2),1),'ub',ones(size(H2_cons.G,2),1));
X_H2_cons = plus(H2_cons.c,affineMap(Box_H2_cons,H2_cons.G));

% Plotting
fig = figure('Position',[100 100 900 600]); hold on
subplot(1,2,1);hold on
hull_uncons.plot('k',0.5);
x.plot('r',1);
y.plot('b',1);
plot(x.c(1),x.c(2),'.k','markersize',20)
plot(y.c(1),y.c(2),'.k','markersize',20)
plot(hull_uncons.c(1),hull_uncons.c(2),'.k','markersize',20)

xlabel('$z_1$')
ylabel('$z_2$')
xlim([-10 2])
ylim([-6 6])
yticks(linspace(-5,5,3))

leg = legend('$Z_h$','$Z_1$','$Z_2$','location','southwest');
set(leg,'Interpreter','latex');
grid off
box on
axis square

subplot(1,2,2);hold on
hull_cons.plot('k',0.5);
x_cons.plot('r',1);
y_cons.plot('b',1);
x.plot('r',0.1);
y.plot('b',0.1);

plot(X_H1_cons,'edgecolor',[.5 0 .5],'linewidth',2)
plot(X_H2_cons,'edgecolor',[.5 0 .5],'linewidth',2)
plot(x_cons.c(1),x_cons.c(2),'.k','markersize',20)
plot(y_cons.c(1),y_cons.c(2),'.k','markersize',20)
plot(hull_cons.c(1),hull_cons.c(2),'.k','markersize',20)

xlabel('$z_1$')
ylabel('$z_2$')
xlim([-10 2])
ylim([-6 6])
yticks(linspace(-5,5,3))
leg = legend('$Z_{ch}$','$Z_{c1}$','$Z_{c2}$','location','southwest')
% leg = legend([p1 p2 p3],'$Z_{ch}$','$Z_{c1}$','$Z_{c2}$','location','southwest');
set(leg,'Interpreter','latex');
grid off
box on
axis square

set(gcf, 'Color', 'w');

% Check 

containCheck(x_cons,hull_cons)
containCheck(y_cons,hull_cons)

%% 7- Computing Hausdorff distance between sets $x_1$ and $x_2$

x1 = conZono;
x1.c = [0;0];
x1.G = 1*eye(2);
x2 = conZono;
x2.c = [0;0];
x2.G = [1 1; 1 -1];

figure; hold on
plot(x1,'r',0.5)
plot(x2,'b',0.5)

xlabel('$z_1$','interpreter','Latex');
ylabel('$z_2$','interpreter','Latex');
set(gcf,'color','w');
grid on;
box on;
legend('$\mathcal{X}_1$','$\mathcal{X}_2$','$\mathcal{X} \cup \mathcal{Y}$','interpreter','Latex','Location','northoutside','NumColumns',2);

d = hausdorffDistance(x1,x2)

%% 8- Constructing and plotting scaled constrained zonotope

figure; hold on
phi = diag([2.5,1,1]);
zAH = conZono2AHPoly(z);
z.G = z.G*phi;
z.plot
zAH.plot('b',0.2)

xlabel('$z_1$','interpreter','Latex');
ylabel('$z_2$','interpreter','Latex');
set(gcf,'color','w');
grid on;
box on;
legend('$\tilde{\mathcal{Z}}$','$\mathcal{Z}$','interpreter','Latex','Location','northoutside','NumColumns',2);

%% 9- Computing the bounds of \xi variables of the constrained zonotope 

x2 = conZono;
x2.c = [0;0];
% x2.G = 1*eye(2);
x2.G = [1 0 1; 1 2 -1];
x2.A = [-2 1 -1];
x2.b = [2];

maxIter = 100
[R,E,rep] = bounds(x2,maxIter)

%% 10- Computing the \xi and \rho bounds for each of the generators considering each constraint separately

x2 = conZono;
x2.c = [0;0];
x2.G = [1 0 1 1; 1 2 -1 3];
x2.A = [-2 1 -1 1; 1 1 -1 1];
x2.b = [2;1];

maxIter = 100
[R,E,rep] = bounds_ind(x2,maxIter) % Need to check

%% 11- Constructing and plotting inner-approximation of a constrained zonotope by a constrained zonotope, zonotope and box

x0 = conZono;
x0.c = [0;0];
x0.G = [-1 3 4;4 -2 -5];
x0.A = zeros(0,3);
x0.b = [];

E = [1 0; 0 1];
f = [5; 5];
H = Polyhedron('H',[E f]);

haus = 0;

[x] = halfspaceIntersection(x0,H); % Halfspace intersection

% 
n_g = size(x.A,2);
n_c = size(x.A,1);

bound_iter = 100;

[R,E,~] = bounds(x,bound_iter);
R_max = max(abs(R)')';
[r_min,j_min] = min(R_max);
[M,i_min] = max(abs(x.A(:,j_min)));

[x_r] = removerowicolumnj(x,i_min,j_min);

x_s = innerApprox(x_r,x,'Infinity Norm');

% Unconstrained zonotope

xpar_T = null(x.A);
xpar_s = pinv(x.A)*x.b;

z = conZono;
z.c = x.c+x.G*xpar_s;
z.G = x.G*xpar_T;
z.A = zeros(0,size(z.G,2));
z.b = [];

z_s = innerApprox(z,x,'Infinity Norm');

% Box 
b = conZono;

b.c = x.c;
b.G = eye(2);
b.A = zeros(0,size(b.G,2));
b.b = [];

b_s = innerApprox(b,x,'Infinity Norm');

CG_ratio = volumeratio(x_s,x);
G_ratio = volumeratio(z_s,z);
B_ratio = volumeratio(b_s,x);

[CG_ratio G_ratio B_ratio]

% Plotting 

figure('Position',[100 100 900 600]); 
subplot(1,2,1);hold on
x.plot('r',1);
x_s.plot('b',1);
% plot(Xs_set,'color','b');

xlim([-8 8])
ylim([-10 6])
xlabel('$z_1$')
ylabel('$z_2$')
yticks(linspace(-10,5,4))
leg = legend('$Z_c$','$Z_r$');
set(leg,'Interpreter','latex','location','northeast');
grid off
box on
axis square

subplot(1,2,2);hold on
x.plot('r',1);
z_s.plot('c',1);
b_s.plot('m',1);

xlim([-8 8])
ylim([-10 6])
xlabel('$z_1$')
ylabel('$z_2$')
yticks(linspace(-10,5,4))
leg = legend('$Z_c$','$Z$','$B$');
set(leg,'Interpreter','latex','location','northeast');
grid off
box on
axis square

%% 12- Computing the Pontryagin difference between a zonotope/constrained zonotope and a zonotope (Iterative method)

n_dof = 3;

Z1 = conZono;
Z1.c = zeros(n_dof,1);
Z1.G = [1 1 0 0; 1 0 1 0; 1 0 0 1];
Z1.A = zeros(0,size(Z1.G,2));
Z1.b = [];

Z2 = conZono;
Z2.c = zeros(n_dof,1);
Z2.G = (1/3)*[-1 1 0 0; 1 0 1 0; 1 0 0 1];
Z2.A = zeros(0,size(Z2.G,2));
Z2.b = [];

% Iterative method
Z_d = pontryagin_diff(Z1,Z2); 

% Approximating Pontryagin difference (1-step method) 
tic
Z_d_approx = pontryagin_diff_1step(Z1,Z2)
toc

figure;
Z_d.plot('b',0.4);
hold on;
Z_d_approx.plot('r',0.2);
xlabel('$z_1$','interpreter','Latex');
ylabel('$z_2$','interpreter','Latex');
set(gcf,'color','w');
grid off;
box on;
legend('$\mathcal{Z}_d$','$\tilde{\mathcal{Z}}_d$','interpreter','Latex','Location','northoutside','NumColumns',2);

volumeratio(Z_d_approx,Z_d)

%% 13- Determining and eliminating redundant generators and constraints of a constrained zonotope

x1 = conZono;
x1.c = [0;0];
x1.G = [1 1;1 -1];
x1.A = [];
x1.b = [];

x2 = conZono;
x2.c = [0;0];
x2.G = [1 0;0 1];
x2.A = [];
x2.b = [];

% Redundant representation
[x] = generalizedIntersection(x1,x2,eye(2));

% Determining redundant generators and constraints
[x_r] = rrefcp(x);
nc = size(x_r.A,1);
ng = size(x_r.A,2);
red_iter = 100;
[redund] = redundancy_indices(x_r,red_iter);

% Removing redundant generators and corresponding constraints
while length(redund) > 0
    [x_r] = removerowicolumnj(x_r,redund(1,1),redund(1,2));
    nc = size(x_r.A,1);
    ng = size(x_r.A,2);
    redund = [];
    [redund] = redundancy_indices(x_r,red_iter);
end

x_r

% Plotting

figure('Position',[100 100 400 600]); hold on
hold on;
x1.plot('r',1);
x_r.plot('b',1);

xlabel('$z_1$')
ylabel('$z_2$')

leg = legend('$Z_1$','$Z_2$');
set(leg,'Interpreter','latex');
grid off
box on
axis square

set(gcf, 'Color', 'w');

%% 14- Reduced-row echelon form (rref) of a constrained zonotope with [A b] in rref.

[out] = rrefcp(Z_d)

%% 15- Rescaling the $\xi$ bounds of a constrained zonotope

[Z_d_rescaled] = rescale(Z_d)

%% 16 - Computing the interval hull of a constrained zonotope

%Define constrained zonotope
x0 = conZono([0;0],[-1 3 4;4 -2 -5],zeros(0,3),[]);
E = [1 0; 0 1];
f = [5; 5];
H = Polyhedron('H',[E f]);
x = halfspaceIntersection(x0,H); % Halfspace intersection

%compute interval hull
[IH,IH_bounds] = x.intervalHull;


figure;
x.plot('b',0.1);hold on;
IH.plot('g',0.1);
xlim([-9 6]);
ylim([-10 6]);
legend('$x$','$IH$','location','southwest');

%% 17 - computing the cartesian product of two constrained zonotopes

z1 = conZono([10;10],[1;1],[3 4; -5 2],[0;0]);
z2 = conZono([0],[5]);

z3 = cartProd(z1,z2)