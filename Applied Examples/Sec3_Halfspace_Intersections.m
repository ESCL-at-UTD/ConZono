%% Plot Settings
set(0,'defaultLineLineWidth', 2)
set(0,'defaultAxesFontName' , 'Times')
set(0,'defaultTextFontName' , 'Times')
set(0,'defaultAxesFontSize' , 18)
set(0,'defaultTextFontSize' , 18)
set(0,'defaulttextinterpreter','latex')
set(0,'defaultAxesGridLineStyle','-.')

%% Section 3, Halfspace Intersections (Fig. 1)
% Generate Zonotope
x.c = [0;0];
x.G = [1 1; 0 2];
x.A = zeros(0,size(x.G,2));
x.b = [];

x.Box = Polyhedron('lb',-ones(size(x.G,2),1),'ub',ones(size(x.G,2),1)); 
x.poly = plus(x.c,affineMap(x.Box,x.G));

%% Generate Halfspace and Hyperplane
E = [3 1];
f = [3];
H_ = Polyhedron('H',[E f]);
H = Polyhedron('He',[E f]);

%% Evaluate intersection
abs(f-E*x.c)
sum(abs(E*x.G))

%% Compute Intersection of Zonotope and Halfspace
 
[z] = halfspaceIntersection(x,H_);

z.Box = Polyhedron('lb',-ones(size(z.G,2),1),'ub',ones(size(z.G,2),1),'He',[z.A z.b]); 
z.poly = plus(z.c,affineMap(z.Box,z.G));

%% Generate Second Halfspace and Hyperplane
E2 = [4 1];
f2 = [6];
H2_ = Polyhedron('H',[E2 f2]);
H2 = Polyhedron('He',[E2 f2]);

%%  Plot Figure

H_mod = Polyhedron('H',[E f],'lb',[-3 -3],'ub',[3 3]);
H2_mod = Polyhedron('H',-[E2 2*f2],'lb',[-3 -3],'ub',[3 3]);

figure('Position',[100 100 900 600])
subplot(1,2,1); hold on
p1 = plot(H_mod,'color',[0.9 0.7 1],'LineStyle','none');
p2 = plot(x.poly,'color',[0.5 0.5 0.5]);
p3 = plot([2 0],[-3 3],'color',[0.6 0.3 1],'LineWidth',3);
p4 = plot(z.poly,'color',[0.4 0.8 0.9]);
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

leg = legend([p2 p4 p1 p3],'$Z$','$Z_h$','$H_-$','$H$');
set(leg,'Interpreter','latex','location','northwest');

subplot(1,2,2); hold on
p1 = plot(H2_mod,'color',[0.9 0.7 1],'LineStyle','none');
p2 = plot(x.poly,'color',[0.5 0.5 0.5]);
p3 = plot([2.25 0.75],[-3 3],'color',[0.6 0.3 1],'LineWidth',3);
p4 = plot(z.poly,'color',[0.4 0.8 0.9]);
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

leg = legend([p2 p4 p1 p3],'$Z$','$Z_c$','$H_-$','$H$');
set(leg,'Interpreter','latex','location','northwest');

% set(gcf, 'Color', 'w');
% export_fig Halfspace_Int.pdf -painters 

%% Check constrained zonotope hyperplane intersection
H2_plus = Polyhedron('H',[-E2 -f2]);

% z_rref = CG_rref(z);
z_rref = z;

% [z_plus] = halfspaceIntersection(z,H2_plus);
[z_plus] = halfspaceIntersection(z_rref,H2_plus);

[z_plus_irred,R_plus,E_plus,~] = Bounds_Rescale(z_plus,2);
% Requires two iterations to check for non-intersection

z_plus_irred
R_plus
E_plus

z_plus_irred_feas = 0;
for i = 1:length(E_plus(:,1))
    if (E_plus(i,1) <= E_plus(i,2))
        z_plus_irred_feas = z_plus_irred_feas + 1;
    end
    if (i == length(E_plus(:,1))&& z_plus_irred_feas == length(E_plus(:,1)))
        disp('Non-empty set');
    elseif (i == length(E_plus(:,1)))
        disp('Empty set');
    end
end
