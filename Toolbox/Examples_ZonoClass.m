%% 1- Constructing and plotting a zonotope 
x = Zono;
x.c = [2;0];
x.G = [1 0;0 1];

c_zon = Zono;
c_zon.c = [2;0];
c_zon.G = [3 0; 0 1];

containCheck(x,c_zon) % Check whether $x \in c$.

figure;
hold on;
c_zon.plot('r',0.1);
x.plot('b',0.2);
xlabel('$z_1$')
ylabel('$z_2$')
grid off
box on
axis square
legend('$\mathcal{X}$','$\mathcal{Z}$','interpreter','Latex')
set(gcf,'color','w');
grid off
box on
axis square

%% 2- Computing and plotting zonotopic inner-approximation 

z = Zono;
z.c = zeros(2,1);
rng(1)
% z.G = rand(2,5);
z.G = [4 3 -2 0.2 0.5; 0 2 3 0.6 -0.3];

z_inapprox = innerApprox(z,3)

figure;
z.plot('r',1)
hold on;
z_inapprox.plot('b',1)
% z.volume
% z_inapprox.volume
volumeratio(z_inapprox,z)
xlabel('$z_1$')
ylabel('$z_2$')
xlim([-10 10])
ylim([-10 10])

grid off
box on
axis square
legend('$\mathcal{Z}$','$\mathcal{Z}_r$','interpreter','Latex');
set(gcf,'color','w')
%% 3- Computing and plotting the minkowski addition of two zonotopes
z1 = Zono;
z1.c = [0;0];
z1.G = [2 1; 1 1];

z2 = Zono;
z2.c = [1;0];
z2.G = [1 2 3; 1 -1 2];

z_plus = plus(z1,z2)
volumeratio(z_inapprox,z)

figure;
hold on;
z1.plot('b',1);
z2.plot('r',0.2);
z_plus.plot('k',0.2);
xlabel('$z_1$')
ylabel('$z_2$')
grid off
box on
axis square
legend('$\mathcal{Z}_1$','$\mathcal{Z}_2$','$\mathcal{Z}_1 \oplus \mathcal{Z}_2$','interpreter','Latex');
set(gcf,'color','w');

%% 4- Scaling and plotting the diagonal elements of a zonotope
z = Zono;
z.c = [1;0];
z.G = [1 2 3; 1 -1 2];

phi = diag([2.5,1,1]);
zAH = Zono2AHPoly(z);
z.G = z.G*phi;
figure;
z.plot
hold on;
zAH.plot('b',0.2)
xlabel('$z_1$')
ylabel('$z_2$')
legend('$\mathcal{Z}$','$\hat{\mathcal{Z}}$','$\mathcal{Z}_1$ \oplus \mathcal{Z}_2','interpreter','Latex')
grid off
box on
axis square
set(gcf,'color','w')
%% 5- Detecting and removing parallel generators and plotting
z = Zono;
z.c = zeros(2,1);
z.G = [1 2 3; 1 2 1];

z_irred  = redundancy_parallelgen(z)

figure;
z.plot
hold on;
z_irred.plot('b',0.2)
xlabel('$z_1$')
ylabel('$z_2$')
legend('$\mathcal{Z}$','$\mathcal{Z}_r$','$\mathcal{Z}_1$ \oplus \mathcal{Z}_2','interpreter','Latex')
grid off
box on
axis square
set(gcf,'color','w')
