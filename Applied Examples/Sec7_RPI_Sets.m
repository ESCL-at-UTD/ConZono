%% Plotting 
set(0,'defaultLineLineWidth', 2)
set(0,'defaultAxesFontName' , 'Times')
set(0,'defaultTextFontName' , 'Times')
set(0,'defaultAxesFontSize' , 18)
set(0,'defaultTextFontSize' , 18)
set(0,'defaulttextinterpreter','latex')
set(0,'defaultAxesGridLineStyle','-.')

%% Section 7, Robust Positively Invariant Sets (Fig. 7)
% Define System
A = [1 1; 0 1];
B = [0.5;1];
n = size(A,1);
m = size(B,2);
Q = eye(2);
W = Polyhedron('lb',-0.1*ones(2,1),'ub',0.1*ones(2,1));

ch = 1; % Select controller case

if (ch == 1) % Case 1: 
    R = 1;
elseif (ch == 2) % Case 2 
    R = 100;
end
[lq_sol,lqr_pol,K] = dare(A,B,Q,R);
K = -K;
Ak = A + B*K; 

%% eps-mRPI Set Outer Approximation (Rakovic et.al)
N_ave = 10;
N_calcs = N_ave+1;

% Using H-Rep
eps_all = logspace(0,-5,6);

mean_tcalc_H_eps_mRPI = zeros(length(eps_all),1);
volume_H_eps_mRPI = zeros(length(eps_all),1);
complexity_H_eps_mRPI = zeros(length(eps_all),1);

% Using G-Rep
mean_tcalc_G_eps_mRPI = zeros(length(eps_all),1);
volume_G_eps_mRPI = zeros(length(eps_all),1);
complexity_G_eps_mRPI = zeros(length(eps_all),1);

for iter = 1:length(eps_all)
    eps = eps_all(iter);
    
    % Using H-Rep
    for i = 1:N_calcs
        tstart = tic;
        H_mRPI_eps = Approx_RPI(Ak,W,eps);
        tcalc_H_eps_mRPI(i) = toc(tstart);
    end
    mean_tcalc_H_eps_mRPI(iter) = mean(tcalc_H_eps_mRPI(2:end));
    volume_H_eps_mRPI(iter) = H_mRPI_eps.volume;
    complexity_H_eps_mRPI(iter) = size(H_mRPI_eps.H,1);
    
    % Using G-Rep
    W_set.c = zeros(2,1);
    W_set.G = 0.1*eye(2);
    
    for i = 1:N_calcs
        tstart = tic;
        G_mRPI_eps = Approx_RPI(Ak,W_set,eps);
        tcalc_G_eps_mRPI(i) = toc(tstart);
    end
    G_mRPI_eps_set = zonotope([G_mRPI_eps.c G_mRPI_eps.G]);
    mean_tcalc_G_eps_mRPI(iter) = mean(tcalc_G_eps_mRPI(2:end));
    volume_G_eps_mRPI(iter) = G_mRPI_eps_set.volume;
    complexity_G_eps_mRPI(iter) = size(G_mRPI_eps.G,2);
end

%% eps_mRPI (High precision)
eps = 10^-9;
G_mRPI_eps = Approx_RPI(Ak,W_set,eps);
G_mRPI_eps_set = zonotope([G_mRPI_eps.c G_mRPI_eps.G]);
volume_best = G_mRPI_eps_set.volume;

%% One Step Method H-Rep ( Trodden )

r = complexity_H_eps_mRPI;

mean_tcalc_H_1_step = zeros(length(r),1);
volume_H_1_step = zeros(length(r),1);
complexity_H_1_step = zeros(length(r),1);

F = W.H(:,1:end-1);
g = W.H(:,end);
for iter = 1:length(r)
    % Compute P s.t. Px <= q
    P = [sin(2*pi*[0:r(iter)-1]'/r(iter)) cos(2*pi*[0:r(iter)-1]'/r(iter))];
    
    % Solving using Yalmip
    yalmip('clear') 
    xi_ = sdpvar(n,r(iter));
    w_ = sdpvar(n,r(iter));
    c_ = sdpvar(r(iter),1);
    d_ = sdpvar(r(iter),1);
    
    objs = 0;
    cons = [];
    for i = 1:r(iter)
        objs = objs - (c_(i) + d_(i));
        cons = [cons, c_(i) <= P(i,:)*Ak*xi_(:,i)];
        cons = [cons, P*xi_(:,i) <= c_ + d_];
        cons = [cons, d_(i) <= P(i,:)*w_(:,i)];
        cons = [cons, F*w_(:,i) <= g];
    end
    opts = sdpsettings('solver','gurobi');
    Opt_solve = optimizer(cons,objs,opts,[],{c_,d_}); % Create exported optimization model object
    
    for i = 1:N_calcs
        tstart = tic;
        [OUT,diagnostics] = Opt_solve([]);
        tcalc_H_1_step(i) = toc(tstart);
    end
    
    c_val = cell2mat(OUT(1));
    d_val = cell2mat(OUT(2));
    q_opt(iter) = mat2cell(c_val + d_val,r(iter),1);
        
    mRPI_ = Polyhedron('H',[P,q_opt{iter}]);

    mean_tcalc_H_1_step(iter) = mean(tcalc_H_1_step(2:end));
    volume_H_1_step(iter) = mRPI_.volume;
    complexity_H_1_step(iter) = r(iter);
    
end

%% One Step Method G-Rep ( New )

n_steps = (complexity_G_eps_mRPI/size(W_set.G,2))-1;

mean_tcalc_G_1_step = zeros(length(n_steps),1);
volume_G_1_step = zeros(length(n_steps),1);
complexity_G_1_step = zeros(length(n_steps),1);

for iter = 1:length(n_steps)
    R = W_set;
    for i = 1:n_steps(iter)
        R.c = R.c + Ak^i*W_set.c;
        R.G = [R.G Ak^i*W_set.G];
    end

    % Optimization for minimum 
    yalmip('clear') 
    c_R = sdpvar(2,1); % Center of mRPI Set 
    S_R = sdpvar(size(R.G,2),1);
    Gamma_s1 = sdpvar(size(R.G,2),size(R.G,2),'full');
    Gamma_s2 = sdpvar(size(R.G,2),size(W_set.G,2),'full');
    beta_s = sdpvar(size(R.G,2),1);

    cons = [];
    cons = [cons, Ak*R.G*diag(S_R) == R.G*Gamma_s1];
    cons = [cons, W_set.G == R.G*Gamma_s2];
    cons = [cons, (eye(2)-Ak)*c_R-W_set.c == R.G*beta_s];
    cons = [cons, abs(Gamma_s1)*ones(size(R.G,2),1)+abs(Gamma_s2)*ones(size(W_set.G,2),1)+abs(beta_s) <= S_R];
    cons = [cons, S_R >= 0];
    objs = norm(S_R,'inf');
    opts = sdpsettings('solver','gurobi');
    Opt_solve = optimizer(cons,objs,opts,[],{c_R,S_R}); % Create exported optimization model object

    for i = 1:N_calcs
        tstart = tic;
        [OUT,diagnostics] = Opt_solve([]);
        tcalc_G_1_step(i) = toc(tstart);
    end
    
    Rs.c = cell2mat(OUT(1));
    Rs.G = R.G*diag(cell2mat(OUT(2)));
    
    Rs_set = zonotope([Rs.c Rs.G]);
    
    mean_tcalc_G_1_step(iter) = mean(tcalc_G_1_step(2:end));
    volume_G_1_step(iter) = Rs_set.volume;
    complexity_G_1_step(iter) = size(Rs.G,2);
    
end
%% Plot eps-mRPI Data
% Plotting H complexity divided by 2 
% Plotting normalized volume error 

volume_ratio_H_eps_mRPI = (volume_H_eps_mRPI./volume_best).^(1/2);
volume_ratio_G_eps_mRPI = (volume_G_eps_mRPI./volume_best).^(1/2);
volume_ratio_H_1_step = (volume_H_1_step./volume_best).^(1/2);
volume_ratio_G_1_step = (volume_G_1_step./volume_best).^(1/2);

%% Plotting
figure('Position',[100 100 900 600]);
f1 = subplot(2,1,1);hold on
f1_pos=get(f1,'Position'); 
set(f1,'Position',[f1_pos(1) f1_pos(2)-0.01 f1_pos(3) f1_pos(4)+0.01])
plot(complexity_H_eps_mRPI/2,volume_ratio_H_eps_mRPI-1,'ro','MarkerSize',10)
plot(complexity_G_eps_mRPI,volume_ratio_G_eps_mRPI-1,'bs','MarkerSize',10)
plot(complexity_H_1_step/2,volume_ratio_H_1_step-1,'k+','MarkerSize',10)
plot(complexity_G_1_step,volume_ratio_G_1_step-1,'mx','MarkerSize',10)
set(f1,'yscale','log')
% title('$\epsilon$-mRPI')
xlabel('$n_g$,$\frac{1}{2}n_h$')
ylabel('$V_r - 1$')
xlim([0 40])
ylim([10^-6 10^1])
yticks([10^-6 10^-3 10^0])
% set(gca,'xticklabel',{[]})
legend('\epsilon-mRPI (H-Rep)','\epsilon-mRPI (G-Rep)','1-step (H-Rep)','1-step (G-Rep)','location','northoutside','orientation','horizontal')
grid on
box on

f2 = subplot(2,1,2);hold on
f2_pos=get(f2,'Position'); 
% set(f2,'Position',[f2_pos(1) f2_pos(2) f2_pos(3)+0.025 f2_pos(4)+0.025])
plot(complexity_H_eps_mRPI/2,mean_tcalc_H_eps_mRPI,'ro','MarkerSize',10)
plot(complexity_G_eps_mRPI,mean_tcalc_G_eps_mRPI,'bs','MarkerSize',10)
plot(complexity_H_1_step/2,mean_tcalc_H_1_step,'k+','MarkerSize',10)
plot(complexity_G_1_step,mean_tcalc_G_1_step,'mx','MarkerSize',10)
set(f2,'yscale','log')
set(gca,'yminorgrid','off')
% title('$\epsilon$-mRPI')
xlabel('$n_g$,$\frac{1}{2}n_h$')
ylabel('$\Delta t_{calc}$ [s]')
xlim([0 40])
ylim([10^-3 10^0])
yticks([10^-3 10^-2 10^-1 10^0])
grid on
box on

set(gcf, 'Color', 'w');
% export_fig RPI_Sets_K1_ave10.pdf -painters
