%% Plot Settings
set(0,'defaultLineLineWidth', 2)
set(0,'defaultAxesFontName' , 'Times')
set(0,'defaultTextFontName' , 'Times')
set(0,'defaultAxesFontSize' , 18)
set(0,'defaultTextFontSize' , 18)
set(0,'defaulttextinterpreter','latex')
set(0,'defaultAxesGridLineStyle','-.')

%% Section 7, Robust Positively Invariant Sets (Fig. 8)
% Scalability with System order for an integrator system 

Sys_ord_array = [2:10];
eps = 10^(-4);
N_ave = 10;
N_calcs = N_ave+1;

ch = 1; % Select controller case

mean_tcalc_H_eps_mRPI = zeros(length(Sys_ord_array),1);
complexity_H_eps_mRPI = zeros(length(Sys_ord_array),1);
mean_tcalc_G_eps_mRPI = zeros(length(Sys_ord_array),1);
complexity_G_eps_mRPI = zeros(length(Sys_ord_array),1);
mean_tcalc_H_1_step = zeros(length(Sys_ord_array),1);
complexity_H_1_step = zeros(length(Sys_ord_array),1);
mean_tcalc_G_1_step = zeros(length(Sys_ord_array),1);
complexity_G_1_step = zeros(length(Sys_ord_array),1); 

for p_ord = 1:size(Sys_ord_array,2)
    ord = Sys_ord_array(p_ord);
    ord
    den = zeros(1,ord+1);
    den(1) = 1;
    sys_c = tf(1,den);
    sys_c = ss(sys_c);
    sys_d = c2d(sys_c,1);  % Convert to discrete time with 1 sec timestep
    A = sys_d.A;
    B = sys_d.B;
    
    n = size(A,1);
    m = size(B,2);
    Q = eye(ord);
    W = Polyhedron('lb',-0.1*ones(ord,1),'ub',0.1*ones(ord,1));
    W_set.c = zeros(ord,1);
    W_set.G = 0.1*eye(ord);
    
    if (ch == 1) % Case 1:
        R = 1;
    elseif (ch == 2) % Case 2
        R = 100;
    end
    [lq_sol,lqr_pol,K] = dare(A,B,Q,R);
    K = -K;
    Ak = A + B*K;

    % eps-mRPI Set Outer Approximation (Rakovic et.al)
%     % Using H-Rep
%     for i = 1:N_calcs
%         tstart = tic;
%         H_mRPI_eps = Approx_RPI(Ak,W,eps);
%         tcalc_H_eps_mRPI(i) = toc(tstart);
%     end
%     mean_tcalc_H_eps_mRPI(p_ord) = mean(tcalc_H_eps_mRPI(2:end));
%     complexity_H_eps_mRPI(p_ord) = size(H_mRPI_eps.H,1);
    
    % Using G-Rep
    for i = 1:N_calcs
        tstart = tic;
        G_mRPI_eps = Approx_RPI(Ak,W_set,eps);
        tcalc_G_eps_mRPI(i) = toc(tstart);
    end
    mean_tcalc_G_eps_mRPI(p_ord) = mean(tcalc_G_eps_mRPI(2:end));
    complexity_G_eps_mRPI(p_ord) = size(G_mRPI_eps.G,2);

    % One step Optimization Technique, RPI set complexity and volume ratio with system order 
    % H-Rep
    if ord <= 6
        %     n_steps = ord-1;
        n_steps = 1;
        R = W;
        for i = 1:n_steps
            i
            R = plus(R,affineMap(W,Ak^i));
            R = minHRep(R);
            size(R.H)
        end
        %     i = 0;
        %     while size(R.H,1) < complexity_H_eps_mRPI(p_ord)
        %         i = i + 1;
        %         R = plus(R,affineMap(W,Ak^i));
        %         R = minHRep(R);
        %     end
        
        r = size(R.H,1);
        F = W.H(:,1:end-1);
        g = W.H(:,end);
        P = R.H(:,1:end-1);
        
        % Solving using Yalmip
        yalmip('clear')
        xi_ = sdpvar(n,r);
        w_ = sdpvar(n,r);
        c_ = sdpvar(r,1);
        d_ = sdpvar(r,1);
        
        objs = 0;
        cons = [];
        for i = 1:r
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
        mean_tcalc_H_1_step(p_ord) = mean(tcalc_H_1_step(2:end));
        complexity_H_1_step(p_ord) = r;
    end
    
    % G-Rep
%     n_steps = (complexity_G_eps_mRPI(p_ord)/ord)-1;
    n_steps = ord-1;
    R = W_set;
    for i = 1:n_steps
        R.c = R.c + Ak^i*W_set.c;
        R.G = [R.G Ak^i*W_set.G];
    end
    
    % Optimization for minimum 
    yalmip('clear') 
    c_R = sdpvar(ord,1); % Center of mRPI Set 
    S_R = sdpvar(size(R.G,2),1);
    Gamma_s1 = sdpvar(size(R.G,2),size(R.G,2),'full');
    Gamma_s2 = sdpvar(size(R.G,2),size(W_set.G,2),'full');
    beta_s = sdpvar(size(R.G,2),1);

    cons = [];
    cons = [cons, Ak*R.G*diag(S_R) == R.G*Gamma_s1];
    cons = [cons, W_set.G == R.G*Gamma_s2];
    cons = [cons, (eye(ord)-Ak)*c_R-W_set.c == R.G*beta_s];
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
    mean_tcalc_G_1_step(p_ord) = mean(tcalc_G_1_step(2:end));
    complexity_G_1_step(p_ord) = size(R.G,2);
    
end
%% Plotting
figure('Position',[100 100 900 600]);
f1 = subplot(2,1,1);hold on
f1_pos=get(f1,'Position'); 
set(f1,'Position',[f1_pos(1) f1_pos(2)-0.01 f1_pos(3) f1_pos(4)+0.01])
% plot(Sys_ord_array,complexity_H_eps_mRPI,'ro','MarkerSize',10)
plot(Sys_ord_array,complexity_G_eps_mRPI,'bs','MarkerSize',10)
plot(Sys_ord_array,complexity_H_1_step,'k+','MarkerSize',10)
plot(Sys_ord_array,complexity_G_1_step,'mx','MarkerSize',10)
set(f1,'yscale','log')
xlabel('$n$')
ylabel('$n_g$, $\frac{1}{2}n_h$')
set(gca,'yminorgrid','off')
ylim([10^0 2*10^3])
xlim([1.5 10.5])
yticks([10^0 10^1 10^2 10^3])
legend('\epsilon-mRPI (G-Rep)','1-step (H-Rep)','1-step (G-Rep)','location','northoutside','orientation','horizontal')
grid on
box on

f2 = subplot(2,1,2);hold on
f2_pos=get(f2,'Position'); 
% set(f2,'Position',[f2_pos(1) f2_pos(2) f2_pos(3)+0.025 f2_pos(4)+0.025])
plot(Sys_ord_array,mean_tcalc_H_eps_mRPI,'ro','MarkerSize',10)
plot(Sys_ord_array,mean_tcalc_G_eps_mRPI,'bs','MarkerSize',10)
plot(Sys_ord_array,mean_tcalc_H_1_step,'k+','MarkerSize',10)
plot(Sys_ord_array,mean_tcalc_G_1_step,'mx','MarkerSize',10)
set(f2,'yscale','log')
set(gca,'yminorgrid','off')
% title('$\epsilon$-mRPI')
xlabel('$n$')
ylabel('$\Delta t_{calc}$ [s]')
ylim([10^-3 10^3])
yticks([10^-3 10^0 10^3])
xlim([1.5 10.5])
grid on
box on

set(gcf, 'Color', 'w');
% export_fig RPI_Sets_Scalability_K1_ave10.pdf
