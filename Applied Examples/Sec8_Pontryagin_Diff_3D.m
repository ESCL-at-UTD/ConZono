%% Plotting 

set(0,'defaultLineLineWidth', 2)
set(0,'defaultAxesFontName' , 'Times')
set(0,'defaultTextFontName' , 'Times')
set(0,'defaultAxesFontSize' , 18)
set(0,'defaultTextFontSize' , 18)
set(0,'defaulttextinterpreter','latex')
set(0,'defaultAxesGridLineStyle','-.')

%% Section 8 - Pontryagin Difference (Fig. 9)
n_dof = 3;

Z1.c = zeros(n_dof,1);
Z1.G = [1 1 0 0; 1 0 1 0; 1 0 0 1];
Z1.A = zeros(0,size(Z1.G,2));
Z1.b = [];

Z2.c = zeros(n_dof,1);
Z2.G = (1/3)*[-1 1 0 0; 1 0 1 0; 1 0 0 1];
Z2.A = zeros(0,size(Z2.G,2));
Z2.b = [];

red_iter = 100;

% Iterative method
[Zint] = Pontryagin_Diff(Z1,Z2,1,red_iter); % Third input is a flag to remove redundancy

Z_d = Zint

T = null(Z_d.A,'r');
s = pinv(Z_d.A)*Z_d.b;
P = Polyhedron('H',[T ones(size(s,1),1)-s; -T ones(size(s,1),1)+s]);
P = minHRep(P);
Z_d_poly = plus(Z_d.c+Z_d.G*s,affineMap(P,Z_d.G*T));
Z_d_poly = minHRep(Z_d_poly);

%% Zonotope Inner-Approx (One-step G-Rep)

[Opt_solve] = Pontryagin_Diff_1step(Z1,Z2);
% [Opt_solve] = Pontryagin_Diff_1step_AH(Z1,Z2);
[OUT,diagnostics] = Opt_solve([]); % Solves the LP
Zd2.c = cell2mat(OUT(1));
S_val = cell2mat(OUT(2));
Zd2.G = [Z1.G Z2.G]*diag(S_val);
Box = Polyhedron('lb',-ones(size(Zd2.G,2),1),'ub',ones(size(Zd2.G,2),1));
Zd2_set = plus(Zd2.c,affineMap(Box,Zd2.G));

%% Plotting Fig. 9
[ratio] = VolumeRatio(Zd2_set,Z_d_poly) % Computes the volume ratio 

figure('Position',[100 100 700 300]); hold on
subplot(1,2,1)
plot(Z_d_poly,'color','r')
view([-20 8]);
xlabel('$z_1$')
ylabel('$z_2$')
zlabel('$z_3$')
xticks(linspace(-1,1,3))
yticks(linspace(-1,1,3))
zticks(linspace(-1,1,3))
leg = legend('$Z_d$');
set(leg,'Interpreter','latex','location','northwest');
grid on
box on
axis square
subplot(1,2,2)
plot(Zd2_set,'color','b')
view([-20 8]);
xlabel('$z_1$')
ylabel('$z_2$')
zlabel('$z_3$')
xticks(linspace(-1,1,3))
yticks(linspace(-1,1,3))
zticks(linspace(-1,1,3))
leg = legend('$\tilde{Z}_d$');
set(leg,'Interpreter','latex','location','northwest');
grid on
box on
axis square

set(gcf, 'Color', 'w');
% export_fig Pontryagin_Diff_3D.pdf -painters 