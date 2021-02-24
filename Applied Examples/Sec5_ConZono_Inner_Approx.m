%% Plot Settings
set(0,'defaultLineLineWidth', 2)
set(0,'defaultAxesFontName' , 'Times')
set(0,'defaultTextFontName' , 'Times')
set(0,'defaultAxesFontSize' , 18)
set(0,'defaultTextFontSize' , 18)
set(0,'defaulttextinterpreter','latex')
set(0,'defaultAxesGridLineStyle','-.')

%% Section 5.2 Inner Approximations - Constrained Zonotopes (Fig. 4)
x0.c = [0;0];
x0.G = [-1 3 4;4 -2 -5];
x0.A = zeros(0,3);
x0.b = [];

x0.Box = Polyhedron('lb',-ones(size(x0.G,2),1),'ub',ones(size(x0.G,2),1));
X0_set = plus(x0.c,affineMap(x0.Box,x0.G));

% E = [4.5 1; -1 2];
% f = [1; 1];
E = [1 0; 0 1];
f = [5; 5];
H = Polyhedron('H',[E f]);

haus = 0;

[x] = halfspaceIntersection(x0,H); % Halfspace intersection

x.Box = Polyhedron('lb',-ones(size(x.G,2),1),'ub',ones(size(x.G,2),1),'He',[x.A x.b]);
X_set = plus(x.c,affineMap(x.Box,x.G));

%%
n_g = size(x.A,2);
n_c = size(x.A,1);

% [x] = CG_rref(x);

bound_iter = 100;

[R,E,~] = Bounds(x,bound_iter);
R_max = max(abs(R)')';
[r_min,j_min] = min(R_max);
[M,i_min] = max(abs(x.A(:,j_min)));

[x_r] = RemoveRowiColumnj(x,i_min,j_min);

% Computing inner-approximating constrained zonotope

[x_s] = ConZono_Containment_Opt(x_r,x); % Computing inner-approximating constrained zonotope

x_s.Box = Polyhedron('lb',-ones(size(x_s.G,2),1),'ub',ones(size(x_s.G,2),1),'He',[x_s.A x_s.b]);
Xs_set = plus(x_s.c,affineMap(x_s.Box,x_s.G));

%% Unconstrained Zonotope
x.T = null(x.A);
x.s = pinv(x.A)*x.b;

z.c = x.c+x.G*x.s;
z.G = x.G*x.T;
z.A = zeros(0,size(z.G,2));
z.b = [];

% Computing inner-approximating zonotope
if haus == 1
    [z_s] = ConZono_Containment_Hausdorff_Opt(z,x);
else
    [z_s] = ConZono_Containment_Opt(z,x);
end

z_s.Box = Polyhedron('lb',-ones(size(z_s.G,2),1),'ub',ones(size(z_s.G,2),1));
Zs_set = plus(z_s.c,affineMap(z_s.Box,z_s.G));

%% Box
b.c = x.c;
b.G = eye(2);
b.A = zeros(0,size(b.G,2));
b.b = [];

haus = 0;
% Computing inner-approximating box

if haus == 1
    [b_s] = ConZono_Containment_Hausdorff_Opt(b,x);
else
    [b_s] = ConZono_Containment_Opt(b,x);
end

b_s.Box = Polyhedron('lb',-ones(size(b_s.G,2),1),'ub',ones(size(b_s.G,2),1));
Bs_set = plus(b_s.c,affineMap(b_s.Box,b_s.G));

%%

CG_ratio = VolumeRatio(Xs_set,X_set);
G_ratio = VolumeRatio(Zs_set,X_set);
B_ratio = VolumeRatio(Bs_set,X_set);

[CG_ratio G_ratio B_ratio]

%%
figure('Position',[100 100 900 600]); 
subplot(1,2,1);hold on
plot(X_set,'color','r');
plot(Xs_set,'color','b');

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
plot(X_set,'color','r');
plot(Zs_set,'color','c');
plot(Bs_set,'color','m');

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

set(gcf, 'Color', 'w');
% export_fig ConZono_Inner_Approx.pdf -painters 
