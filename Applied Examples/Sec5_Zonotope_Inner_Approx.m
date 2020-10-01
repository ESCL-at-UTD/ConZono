%% Section 5 - Inner Approximations - Example 3(Fig. 3)
% Plot Settings
set(0,'defaultLineLineWidth', 2)
set(0,'defaultAxesFontName' , 'Times')
set(0,'defaultTextFontName' , 'Times')
set(0,'defaultAxesFontSize' , 18)
set(0,'defaultTextFontSize' , 18)
set(0,'defaulttextinterpreter','latex')
set(0,'defaultAxesGridLineStyle','-.')

%% Generating Z

Z.c = [0;0];
Z.G = [4 3 -2 0.2 0.5; 0 2 3 0.6 -0.3];
n_g = size(Z.G,2);
n_r = 3;

Box = Polyhedron('lb',-ones(size(Z.G,2),1),'ub',ones(size(Z.G,2),1)); 
Z_poly = plus(Z.c,affineMap(Box,Z.G));

%% Sorting generators and computing T
% Compute generator norms
Norms = vecnorm(Z.G);
% Sort generator norms in decreasing order
[Norms_sort,indx_sort] = sort(Norms,'descend');
Z.G = Z.G(:,indx_sort);

% Seperate n_r largest generators
Z1.G = Z.G(:,1:n_r);
Z2.G = Z.G(:,n_r+1:end);

% Compute magnitude of dot product between generators in Z1 and Z2
alpha_abs = zeros(n_r,n_g-n_r);
alpha = zeros(n_r,n_g-n_r);
for i = 1:n_r
    for j = 1:n_g-n_r
        alpha_abs(i,j) = abs(Z1.G(:,i)'*Z2.G(:,j));
        alpha(i,j) = Z1.G(:,i)'*Z2.G(:,j);
    end
end
% Normalize the dot product with respect to largest dot product
alpha_norm = alpha*diag(1./max(alpha_abs,[],1));

T2 = zeros(n_r,n_g-n_r);
T2(alpha_norm==1) = 1;
T2(alpha_norm==-1) = -1;

T = [eye(n_r);T2'];

%% Computing reduced zonotope $Z_r$
Zr.c = Z.c;
Zr.G = Z.G*T;

Box = Polyhedron('lb',-ones(size(Zr.G,2),1),'ub',ones(size(Zr.G,2),1)); 
Zr_poly = plus(Zr.c,affineMap(Box,Zr.G));

ratio = VolumeRatio(Zr_poly,Z_poly)

%% Plotting

figure('Position',[100 100 400 600]); hold on
plot(Z_poly,'color','r')
plot(Zr_poly,'color','b')

xlim([-10 10])
ylim([-10 10])
xlabel('$z_1$')
ylabel('$z_2$')
yticks(linspace(-10,10,3))
leg = legend('$Z$','$Z_r$');
set(leg,'Interpreter','latex','location','northeast','Orientation','horizontal');
grid off
box on
axis square

set(gcf, 'Color', 'w');
% export_fig Zono_Inner_Approx_Tmat.pdf -painters 

%% Random Zonotopes - Inner Approximation 

n_rand = 100;

for k= 1:n_rand
    rng(k);
    Z.c = rand(2,1);
    rng(k+1);
    Z.G = rand(2,5);
    n_g = size(Z.G,2);
    n_r = 3;

    Box = Polyhedron('lb',-ones(size(Z.G,2),1),'ub',ones(size(Z.G,2),1)); 
    Z_poly = plus(Z.c,affineMap(Box,Z.G));

    % Compute generator norms
    Norms = vecnorm(Z.G);
    % Sort generator norms in decreasing order
    [Norms_sort,indx_sort] = sort(Norms,'descend');
    Z.G = Z.G(:,indx_sort);

    % Seperate n_r largest generators
    Z1.G = Z.G(:,1:n_r);
    Z2.G = Z.G(:,n_r+1:end);

    % Compute magnitude of dot product between generators in Z1 and Z2
    alpha_abs = zeros(n_r,n_g-n_r);
    alpha = zeros(n_r,n_g-n_r);
    for i = 1:n_r
        for j = 1:n_g-n_r
            alpha_abs(i,j) = abs(Z1.G(:,i)'*Z2.G(:,j));
            alpha(i,j) = Z1.G(:,i)'*Z2.G(:,j);
        end
    end
    % Normalize the dot product with respect to largest dot product
    alpha_norm = alpha*diag(1./max(alpha_abs,[],1));

    T2 = zeros(n_r,n_g-n_r);
    T2(alpha_norm==1) = 1;
    T2(alpha_norm==-1) = -1;

    T = [eye(n_r);T2'];

    Zr.c = Z.c;
    Zr.G = Z.G*T;

    Box = Polyhedron('lb',-ones(size(Zr.G,2),1),'ub',ones(size(Zr.G,2),1)); 
    Zr_poly = plus(Zr.c,affineMap(Box,Zr.G));

    ratio(k) = VolumeRatio(Zr_poly,Z_poly)
end

ratio_avg = mean(ratio)