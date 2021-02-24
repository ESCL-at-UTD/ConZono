%% Plot Settings
set(0,'defaultLineLineWidth', 2)
set(0,'defaultAxesFontName' , 'Times')
set(0,'defaultTextFontName' , 'Times')
set(0,'defaultAxesFontSize' , 18)
set(0,'defaultTextFontSize' , 18)
set(0,'defaulttextinterpreter','latex')
set(0,'defaultAxesGridLineStyle','-.')

%% Section 6 Convex Hull
% Zonotopes

x.c = [0;0];
x.G = [0 1 0; 1 1 2];
x.A = zeros(0,3);
x.b = [];

Box_uncons = Polyhedron('lb',-ones(size(x.G,2),1),'ub',ones(size(x.G,2),1));
X_uncons = plus(x.c,affineMap(Box_uncons,x.G));

y.c = [-5;0];
y.G = 0.5*[-1 2 -4;1 1 3];
y.A = zeros(0,3);
y.b = [];

Box_uncons = Polyhedron('lb',-ones(size(y.G,2),1),'ub',ones(size(y.G,2),1));
Y_uncons = plus(y.c,affineMap(Box_uncons,y.G));

hull_uncons = cvxhull(x,y);
size(hull_uncons.A)

Box = Polyhedron('lb',-ones(size(hull_uncons.G,2),1),'ub',ones(size(hull_uncons.G,2),1),'He',[hull_uncons.A hull_uncons.b]);
Hull_uncons = plus(hull_uncons.c,affineMap(Box,hull_uncons.G));

%% Constrained Zonotopes (Calculating cvxhull)
% Intersecting x with a halfspace 
E1 = [1 1];
f1 = [0];
H1 = Polyhedron('H',[E1 f1]);

[x_cons] = halfspaceIntersection(x,H1);
Box = Polyhedron('lb',-ones(size(x_cons.G,2),1),'ub',ones(size(x_cons.G,2),1),'He',[x_cons.A x_cons.b]);
X_cons = plus(x_cons.c,affineMap(Box,x_cons.G));

% Intersecting y with a halfspace 
E2 = [-2.5 1];
f2 = 9.5;
H2 = Polyhedron('H',[E2 f2]);

[y_cons] = halfspaceIntersection(y,H2);
Box = Polyhedron('lb',-ones(size(y_cons.G,2),1),'ub',ones(size(y_cons.G,2),1),'He',[y_cons.A y_cons.b]);
Y_cons = plus(y_cons.c,affineMap(Box,y_cons.G));

%% Computing the cvxhull
hull_cons = cvxhull(x_cons,y_cons);

Box = Polyhedron('lb',-ones(size(hull_cons.G,2),1),'ub',ones(size(hull_cons.G,2),1),'He',[hull_cons.A hull_cons.b]);
Hull_cons = plus(hull_cons.c,affineMap(Box,hull_cons.G));

% Constraining the length of the halfspaces H1, H2, H3,H4 for easy plotting
H1_cons.G = 1.25*[-1.5 1.5]' ;
H1_cons.c = [0;0];
Box_H1_cons = Polyhedron('lb',-ones(size(H1_cons.G,2),1),'ub',ones(size(H1_cons.G,2),1));
X_H1_cons = plus(H1_cons.c,affineMap(Box_H1_cons,H1_cons.G));
H2_cons.c = [-4 -0.5]' ;
H2_cons.G = 1.25*[0.75 ; 1.25*1.5];
Box_H2_cons = Polyhedron('lb',-ones(size(H2_cons.G,2),1),'ub',ones(size(H2_cons.G,2),1));
X_H2_cons = plus(H2_cons.c,affineMap(Box_H2_cons,H2_cons.G));

%% Plotting
fig = figure('Position',[100 100 900 600]); hold on
subplot(1,2,1);hold on
plot(Hull_uncons,'color',[0.5 0.5 0.5],'linewidth',2)
plot(X_uncons,'color','r')
plot(Y_uncons,'color','b')
plot(Hull_uncons,'alpha',0,'linewidth',3)
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
p1 = plot(Hull_cons,'color',[0.5 0.5 0.5],'linewidth',2);
plot(X_uncons,'color','r','alpha',0.1)
plot(Y_uncons,'color','b','alpha',0.1)
p2 = plot(X_cons,'color','r');
p3 = plot(Y_cons,'color','b');
plot(X_H1_cons,'edgecolor',[.5 0 .5],'linewidth',2)
plot(X_H2_cons,'edgecolor',[.5 0 .5],'linewidth',2)
plot(Hull_cons,'alpha',0,'linewidth',3)
plot(x_cons.c(1),x_cons.c(2),'.k','markersize',20)
plot(y_cons.c(1),y_cons.c(2),'.k','markersize',20)
plot(hull_cons.c(1),hull_cons.c(2),'.k','markersize',20)

xlabel('$z_1$')
ylabel('$z_2$')
xlim([-10 2])
ylim([-6 6])
yticks(linspace(-5,5,3))

leg = legend([p1 p2 p3],'$Z_{ch}$','$Z_{c1}$','$Z_{c2}$','location','southwest');
set(leg,'Interpreter','latex');
grid off
box on
axis square

set(gcf, 'Color', 'w');
% export_fig Convex_Hull.pdf -painters 

