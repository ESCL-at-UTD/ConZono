%% Plot Settings
set(0,'defaultLineLineWidth', 2)
set(0,'defaultAxesFontName' , 'Times')
set(0,'defaultTextFontName' , 'Times')
set(0,'defaultAxesFontSize' , 18)
set(0,'defaultTextFontSize' , 18)
set(0,'defaulttextinterpreter','latex')
set(0,'defaultAxesGridLineStyle','-.')

%% Section 5 - Inner Approximations (Fig. 5), Random examples.
N = 100;
haus = 0;
V_ratios = zeros(N,3);
for iter = 1:N
    
rng(iter)

n_g = randi([4 20],1,1);
n_c = randi([1 floor(n_g/2)],1,1);
    
x0.c = [0;0];
x0.G = rand(2,n_g);
x.c = x0.c;
x.G = x0.G;
x.A = rand(n_c,n_g);
x.b = rand(n_c,1);

x.Box = Polyhedron('lb',-ones(size(x.G,2),1),'ub',ones(size(x.G,2),1),'He',[x.A x.b]);
X_set = plus(x.c,affineMap(x.Box,x.G));

%
n_g = size(x.A,2);
n_c = size(x.A,1);

% [x] = CG_rref(x);
bounds_iter = 100;
[R,E] = Bounds(x,bounds_iter);
R_max = max(abs(R)')';
[r_min,j_min] = min(R_max);
[M,i_min] = max(abs(x.A(:,j_min)));

[x_r] = RemoveRowiColumnj(x,i_min,j_min);

[x_s] = ConZono_Containment_Opt(x_r,x);

x_s.Box = Polyhedron('lb',-ones(size(x_s.G,2),1),'ub',ones(size(x_s.G,2),1),'He',[x_s.A x_s.b]);
Xs_set = plus(x_s.c,affineMap(x_s.Box,x_s.G));

% Unconstrained Zonotope
x.T = null(x.A);
x.s = pinv(x.A)*x.b;

z.c = x.c+x.G*x.s;
z.G = x.G*x.T;
z.A = zeros(0,size(z.G,2));
z.b = [];

if haus == 1
    [z_s] = ConZono_Containment_Hausdorff_Opt(z,x);
else
    [z_s] = ConZono_Containment_Opt(z,x);
end

z_s.Box = Polyhedron('lb',-ones(size(z_s.G,2),1),'ub',ones(size(z_s.G,2),1));
Zs_set = plus(z_s.c,affineMap(z_s.Box,z_s.G));

% Box
b.c = x.c;
b.G = eye(2);
b.A = zeros(0,size(b.G,2));
b.b = [];

if haus == 1
    [b_s] = ConZono_Containment_Hausdorff_Opt(b,x);
else
    [b_s] = ConZono_Containment_Opt(b,x);
end


b_s.Box = Polyhedron('lb',-ones(size(b_s.G,2),1),'ub',ones(size(b_s.G,2),1));
Bs_set = plus(b_s.c,affineMap(b_s.Box,b_s.G));

%
% figure; hold on
% % plot(X0_set,'color','r','alpha',0.2);
% plot(X_set,'color','r');
% plot(Xs_set,'color','b');
% plot(Zs_set,'color','k');
% plot(Bs_set,'color','m');

CG_ratio = VolumeRatio(Xs_set,X_set);
G_ratio = VolumeRatio(Zs_set,X_set);
B_ratio = VolumeRatio(Bs_set,X_set);

V_ratios(iter,:) = [CG_ratio G_ratio B_ratio];
end

figure('Position',[100 100 400 300]); 
boxplot(V_ratios,'Labels',{'$Z_r$','$Z$','$B$'})
set(gca,'TickLabelInterpreter','latex')
ylabel('$V_r$')
ylim([0 1.1])
yticks(linspace(0,1,6))
grid off
box on

set(gcf, 'Color', 'w');

% export_fig ConZono_Inner_Approx_Rand.pdf -painters 
