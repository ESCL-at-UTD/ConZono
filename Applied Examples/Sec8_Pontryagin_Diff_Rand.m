set(0,'defaultLineLineWidth', 2)
set(0,'defaultAxesFontName' , 'Times')
set(0,'defaultTextFontName' , 'Times')
set(0,'defaultAxesFontSize' , 18)
set(0,'defaultTextFontSize' , 18)
set(0,'defaulttextinterpreter','latex')
set(0,'defaultAxesGridLineStyle','-.')
%%
n_runs = 100;

% Set dimension and set complexity as in Table 1.
n_ng = [2 4 4;...
        2 8 4;...
        2 4 8;...
        2 8 8;...
        3 6 6;...
        3 12 6;...
        3 6 12;...
        3 12 12;...
        4 8 8;...
        4 16 8;...
        4 8 16;...
        4 16 16];

%% Final Z_d complexity 

for i = 1:size(n_ng,1)
    n_dim = n_ng(i,1);
    n_c1 = 0;
    n_g1 = n_ng(i,2);
    n_g2 = n_ng(i,3);
    n_gd = ((2)^(n_g2))*n_g1;
    n_cd = ((2)^(n_g2))*n_c1 + ((2)^(n_g2) - 1)*n_dim;
    Z_d_cmplx(i,:) = [n_cd n_gd];
end
disp(Z_d_cmplx)

%% Compute Pontryagin difference
n_cases = size(n_ng,1);

Results_Table = zeros(n_cases,11);


for row = 1:n_cases 
    iter = 1;
    i = 0;
    tcalc_diff_CG = zeros(n_runs,1);
    nc_ng_CG = zeros(n_runs,2);
    tcalc_diff_H = zeros(n_runs,1);
    nH = zeros(n_runs,1);
    tcalc_diff_G = zeros(n_runs,1);
    ng_G = zeros(n_runs,1);
    Vr_G = zeros(n_runs,1);
    while iter <= n_runs
        rng(iter+i)
        Lm = 3*n_ng(row,3)/n_ng(row,2);
        Ls = 1;
        [Zm] = randomZonotope(n_ng(row,1),n_ng(row,2),'uniform',Lm);
        [Zs] = randomZonotope(n_ng(row,1),n_ng(row,3),'uniform',Ls);
        z_m.c = Zm.Z(:,1);
        z_m.G = Zm.Z(:,2:end);
        z_m.A = sparse(0,n_ng(row,2));
        z_m.b = [];
        z_s.c = Zs.Z(:,1);
        z_s.G = Zs.Z(:,2:end);
        z_s.A = sparse(0,n_ng(row,3)); % Use to avoid memory limitations
        z_s.b = [];
        Z1 = z_m;
        Z2 = z_s;
        Box_Z1 = Polyhedron('lb',-ones(size(Z1.G,2),1),'ub',ones(size(Z1.G,2),1));
        Z1_HRep = plus(Z1.c,affineMap(Box_Z1,Z1.G));
        Box_Z2 = Polyhedron('lb',-ones(size(Z2.G,2),1),'ub',ones(size(Z2.G,2),1));
        Z2_HRep = plus(Z2.c,affineMap(Box_Z2,Z2.G));
        
        % H-Rep
        tstart = tic;
        Zd_HRep = minus(Z1_HRep,Z2_HRep);
        tcalc_diff_H(iter) = toc(tstart);
        nH(iter) = size(Zd_HRep.H,1);
        try
            volume_H = Zd_HRep.volume;
        catch
            i = i+1;
            continue
        end
        if volume_H < 0.1
            i = i+1;
            continue
        end
                
        % CG-Rep
        tstart = tic;
        [z_d] = Pontryagin_Diff(z_m,z_s,0); % Third input is a flag to remove redundancy
        tcalc_diff_CG(iter) = toc(tstart);
        nc_ng_CG(iter,:) = size(z_d.A);
        
        % G-Rep
        [Opt_solve] = Pontryagin_Diff_1step(Z1,Z2);   
%         [Opt_solve] = Pontryagin_Diff_1step_AH(Z1,Z2);
        tstart = tic;
        [OUT,diagnostics] = Opt_solve([]);
        Zd2.c = cell2mat(OUT(1));
        S_val = cell2mat(OUT(2));
        Zd2.G = [Z1.G Z2.G]*diag(S_val);
        % Remove all zero generators
        Gen_Norms = vecnorm(Zd2.G); 
        not_Zeros = Gen_Norms >= 1e-3;  %%%%%%%%%%%% Removing generators with norm < 1e-3
        Zd2.G = Zd2.G(:,not_Zeros);
        tcalc_diff_G(iter) = toc(tstart);
        Box_Z = Polyhedron('lb',-ones(size(Zd2.G,2),1),'ub',ones(size(Zd2.G,2),1));
        Zd2_HRep = plus(Zd2.c,affineMap(Box_Z,Zd2.G));
        try
            volume_G = Zd2_HRep.volume;
        catch
            i = i+1;
            continue
        end
%         if volume_G < 0.1
%             pause
%             i = i+1;
%             continue
%         end
        ng_G(iter) = size(Zd2.G,2);    
        Vr_G(iter) = VolumeRatio(Zd2_HRep,Zd_HRep);
        
        [row iter]
        iter = iter + 1;
    end
    % Store results
    Results_Table(row,1:3) = n_ng(row,:);
    Results_Table(row,4) = mean(nH);
    Results_Table(row,5) = mean(tcalc_diff_H);
    Results_Table(row,6:7) = [mean(nc_ng_CG(:,1)) mean(nc_ng_CG(:,2))];
    Results_Table(row,8) = Results_Table(row,5)/mean(tcalc_diff_CG);
    Results_Table(row,9) = mean(ng_G);
    Results_Table(row,10) = mean(Vr_G);
    Results_Table(row,11) = Results_Table(row,5)/mean(tcalc_diff_G);
    
end
