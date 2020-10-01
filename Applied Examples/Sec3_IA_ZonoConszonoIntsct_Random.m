%% Plot Settings
set(0,'defaultLineLineWidth', 2)
set(0,'defaultAxesFontName' , 'Times')
set(0,'defaultTextFontName' , 'Times')
set(0,'defaultAxesFontSize' , 18)
set(0,'defaultTextFontSize' , 18)
set(0,'defaulttextinterpreter','latex')
set(0,'defaultAxesGridLineStyle','-.')

plotFlag = 0;
%% Section 3 - Halfspace Intersection, Remark 2
% Generate Zonotope
x.c = [0;0];
x.G = [1 1; 0 2];
x.A = zeros(0,size(x.G,2));
x.b = [];

if plotFlag == 1
    x.Box = Polyhedron('lb',-ones(size(x.G,2),1),'ub',ones(size(x.G,2),1));
    x.poly = plus(x.c,affineMap(x.Box,x.G));
end
%% Random Intersection 
E1 = [3 1];
f1 = 3;
H_ = Polyhedron('H',[E1 f1]);
[z] = halfspaceIntersection(x,H_); % Constrained zonotope

n_runs = 100;
pass = inf(n_runs,1);  % Output will be infinity if any random halfspaces do not intersect the parent zonotope
intersected = inf(n_runs,1);
for i = 1:n_runs
    rng(i)
    max_slope = 10;
    slope = -max_slope*rand;
    point = [1;1] + rand(2,1);
    E2_rnd = [1 -slope];
    f2_rnd = point(2)-slope*point(1); % Point slope formula to create "random" lines
    H2_ = Polyhedron('H',[E2_rnd f2_rnd]);
    
    [z_2] = halfspaceIntersection(x,H2_); % Intersection test using parent zonotope x
    
    if (size(z_2.A,1) ~= 0)  % Parent zonotope x does intersect
        
%         z_LP = conszonohalfspaceIntersection_LP(z,H2_.H); % Check intersection with LP
        z_LP = conszonohalfspaceIntersection_LP(z,H2_); % Check intersection with LP
                
        z_rref = CG_rref(z);
        
        H2_plus = Polyhedron('H',[-E2_rnd -f2_rnd]);
        
        [z_plus] = halfspaceIntersection(z,H2_plus);
        
        maxIter = 100;
        [R_plus,E_plus,rep(i)] = Bounds(z_plus,maxIter); % Iteratively computes E until convergence or max iter
        
        emptySet = max(E_plus(:,1) > E_plus(:,2)); % One if detects empty set -> no intersection
%         [i emptySet]
        
        if (size(z_LP.A,1) == 1) % Constrained zonotope z does not intersect
            pass(i) = emptySet == 1;
            intersected(i) = 0;
        else                     % Constrained zonotope z does intersect
            pass(i) = emptySet == 0;
            intersected(i) = 1;
        end
        
        if plotFlag == 1 % Plot
            z.Box = Polyhedron('lb',-ones(size(z.G,2),1),'ub',ones(size(z.G,2),1),'He',[z.A z.b]);
            z.poly = plus(z.c,affineMap(z.Box,z.G));
            H2_mod = Polyhedron('H',[E2_rnd f2_rnd],'lb',[-3 -3],'ub',[3 3]);
            figure;
            p1 = plot(H2_mod,'color',[0.9 0.7 1],'LineStyle','none');
            hold on;
            p2 = plot(z.poly,'color',[0.4 0.8 0.9]);
            p3 = plot(x.poly,'color','g','alpha',0.2);
            set(gca, 'Layer','top')
            z_plus.Box = Polyhedron('lb',-ones(size(z_plus.G,2),1),'ub',ones(size(z_plus.G,2),1),'He',[z_plus.A z_plus.b]);
            z_plus.poly = plus(z_plus.c,affineMap(z_plus.Box,z_plus.G));
            p4 = plot(z_plus.poly,'color',[0.4 0.8 0.9]);
            z_plus.Box = Polyhedron('lb',-ones(size(z_plus.G,2),1),'ub',ones(size(z_plus.G,2),1),'He',[z_plus.A z_plus.b]);
            z_plus.poly = plus(z_plus.c,affineMap(z_plus.Box,z_plus.G));
            p5 = plot(z_plus.poly,'color','k');
            drawnow
        end        
    end
end
[n_runs sum(pass) sum(intersected) mean(rep)] 

rep_val_non_1 = []; % Number of instances requiring more than 1 iteration to refine E.
for i = 1:n_runs 
    if (rep(i) ~= 1)
        rep_val_non_1 = [rep_val_non_1;rep(i)]
    end
end

sum(rep_val_non_1)
length(rep_val_non_1)