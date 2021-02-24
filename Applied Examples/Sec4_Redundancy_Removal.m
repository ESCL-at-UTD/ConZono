%% Plot Settings
set(0,'defaultLineLineWidth', 2)
set(0,'defaultAxesFontName' , 'Times')
set(0,'defaultTextFontName' , 'Times')
set(0,'defaultAxesFontSize' , 18)
set(0,'defaultTextFontSize' , 18)
set(0,'defaulttextinterpreter','latex')
set(0,'defaultAxesGridLineStyle','-.')

%% Detecting and adding parallel generators of a Zonotope
Z.c = [0;0]; 
scalar = 2;
Z.G = [scalar*(1+1e-3) 1 1; scalar*2 2 -1];

ng = size(Z.G,2);
eps = 1e-6; % Precision

for i = 1:ng
    for j = i+1:ng
        abs(dot(Z.G(:,i),Z.G(:,j)))/(norm(Z.G(:,i))*norm(Z.G(:,j))) >= 1 - eps
    end
end

%% Detecting and adding parallel generators of a constrained Zonotope
Z.c = [0;0]; 
scalar = 2;
Z.G = [scalar*(1+1e-3) 1 1 1; scalar*2 2 -1 3];
Z.A = [8 4 1 0];
Z.b = [1];

Zlift.G = [Z.G;Z.A]; % Lifted zonotope
Zlift.c = [Z.c;-Z.b];

ng = size(Zlift.G,2);
eps = 1e-6; % Precision

for i = 1:ng
    for j = i+1:ng
        abs(dot(Zlift.G(:,i),Zlift.G(:,j)))/(norm(Zlift.G(:,i))*norm(Zlift.G(:,j))) >= 1 - eps
    end
end

%% Section 4 - Redundancy Removal - Example 2(Fig. 2)
% Create and plot a redundant zonotope
x1.c = [0;0];
x1.G = [1 1;1 -1];
x1.A = [];
x1.b = [];

Box1 = Polyhedron('lb',-ones(size(x1.G,2),1),'ub',ones(size(x1.G,2),1));
X1 = plus(x1.c,affineMap(Box1,x1.G));

x2.c = [0;0];
x2.G = [1 0;0 1];
x2.A = [];
x2.b = [];

Box2 = Polyhedron('lb',-ones(size(x2.G,2),1),'ub',ones(size(x2.G,2),1));
X2 = plus(x2.c,affineMap(Box2,x2.G));

% Redundant 
[x] = generalizedIntersection(x1,x2,eye(2));

Box = Polyhedron('lb',-ones(size(x.G,2),1),'ub',ones(size(x.G,2),1),'He',[x.A x.b]);
X = plus(x.c,affineMap(Box,x.G));

%% Determine redundant generators and constraints
x_r = CG_rref(x);
nc = size(x_r.A,1);
ng = size(x_r.A,2);
red_iter = 100;
[redund] = Redundancy_Indices(x_r,red_iter);

%% Remove redundant generators and corresponding constraints
while length(redund) > 0
    [x_r] = RemoveRowiColumnj(x_r,redund(1,1),redund(1,2));
    nc = size(x_r.A,1);
    ng = size(x_r.A,2);
    redund = [];
    [redund] = Redundancy_Indices(x_r,red_iter);
end

Box_r = Polyhedron('lb',-ones(size(x_r.G,2),1),'ub',ones(size(x_r.G,2),1),'He',[x_r.A x_r.b]);
X_r = plus(x_r.c,affineMap(Box_r,x_r.G));

x_r

%% Plot

figure('Position',[100 100 400 600]); hold on
plot(X1,'color','r','alpha',1)
plot(X2,'color','b')

xlabel('$z_1$')
ylabel('$z_2$')

leg = legend('$Z_1$','$Z_2$');
set(leg,'Interpreter','latex');
grid off
box on
axis square

set(gcf, 'Color', 'w');
% export_fig Redundancy_Removal.pdf -painters 

%% Redundancy Removal (X \subseteq Y) - Random generation of Z_1 and Z_2

% Create and plot a redundant zonotope
n_iter = 100;
size_all = zeros(n_iter,2);
size_red_op = zeros(n_iter,1);
count = 0;

for i = 1:n_iter
    
x1.c = [0;0];
%rng(i+1);
x1.G = [1 1;1 -1];
x1.A = [];
x1.b = [];

Box1 = Polyhedron('lb',-ones(size(x1.G,2),1),'ub',ones(size(x1.G,2),1));
X1 = plus(x1.c,affineMap(Box1,x1.G));

x2.c = [0;0];
rng(i);
x2.G = rand(2,2);
x2.A = [];
x2.b = [];

Box2 = Polyhedron('lb',-ones(size(x2.G,2),1),'ub',ones(size(x2.G,2),1));
X2 = plus(x2.c,affineMap(Box2,x2.G));

if (X1.contains(X2))
    count = count + 1;
    % Redundant 
    [x] = generalizedIntersection(x1,x2,eye(2));

    Box = Polyhedron('lb',-ones(size(x.G,2),1),'ub',ones(size(x.G,2),1),'He',[x.A x.b]);
    X = plus(x.c,affineMap(Box,x.G));

    size(x.A);
    % Determine redundant generators
    x_r = CG_rref(x);
    nc = size(x_r.A,1);
    ng = size(x_r.A,2);
    [redund] = Redundancy_Indices(x_r,red_iter);
%
    while length(redund) > 0
        [x_r] = RemoveRowiColumnj(x_r,redund(1,1),redund(1,2));
        nc = size(x_r.A,1);
        ng = size(x_r.A,2);
        [redund] = Redundancy_Indices(x_r,red_iter);
    end
    size_all(i,:) = size(x_r.A);
    if(size(x_r.A,1) == 0)
       size_red_op(i) = 1;
    else
       size_red_op(i) = 0;
    end        
end
% Box_r = Polyhedron('lb',-ones(size(x_r.G,2),1),'ub',ones(size(x_r.G,2),1),'He',[x_r.A x_r.b]);
% X_r = plus(x_r.c,affineMap(Box_r,x_r.G));

%x_r
end

disp(count)
disp(sum(size_red_op,1))