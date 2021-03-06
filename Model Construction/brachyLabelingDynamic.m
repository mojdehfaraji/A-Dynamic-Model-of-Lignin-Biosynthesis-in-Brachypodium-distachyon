function [ ] = brachyLabelingDynamicSmartSearch()
    clc
    
    load StdStt&P_Results;
    p_l = P_val;

    s = 0;

    r = 0.90;
    rc = r / (1 - r);
    
    [S,S_p,D_xp,D_xn,pools,S_x] = Stoi(r);
    
    load StdSttsplrForDynamic1
    StdStt_l = StdStt;
    [r_StdStt,~,~] = size(StdStt_l);
    w = zeros(1,r_StdStt);
        
    flux = 56;
    metabolite = 34;
    
    CA_U = 130;
    PCA500_U = 30;
    
        %1	2	3	4	5	6	7	8	9	10	11	12	13	14	15	16	17	18	19	20	21	22	23	24	25	26	27	28	29	30	31	32	33	34
    I = [0	0	0	0	0	0	0	0	53	0	47	0	0	0	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   0]';

        %       phe_L  tyr_L   G_L     S_L    tot_L   pca_L    fa_L    H_L
    labeling = [0.35    0.0    0.223   0.210  0.222   0.21     0.23    0.0   % labelled phe CTRL
                0.0     0.35   0.165   0.181  0.186   0.17     0.13    0.0   % labelled tyr CTRL
                1.0     0.0    0.576   0.495  0.534   0.00     0.00    0.65  % all labelled phe
                0.0     1.0    0.424   0.505  0.466   0.00     0.00    0.35];% labelled tyr CTRL

    I_L = I .* [0	0	0	0	0	0	0	0	labeling(1,1)	0	labeling(1,2)	0	0	0	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   0
                0	0	0	0	0	0	0	0	labeling(2,1)	0	labeling(2,2)	0	0	0	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   0
                0	0	0	0	0	0	0	0	labeling(3,1)	0	labeling(3,2)	0	0	0	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   0
                0	0	0	0	0	0	0	0	labeling(4,1)	0	labeling(4,2)	0	0	0	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   0]';
    I_U = I - I_L;
          
    tol = 0.2;

    n_v = 38;
    n_e = 4;
    n_d = 14;
    n_h = 3;
    n_p = n_v + n_e + 4 + n_h;

    [org_sample,c_p_l] = size(p_l);

    admsbl_p = zeros(1,c_p_l,r_StdStt);
    slicer = ones(1,r_StdStt);
    admsbl_pSliced = mat2cell(admsbl_p,1,c_p_l,slicer);
    
    V_results = double.empty;
    X_results = double.empty;
    P_results = double.empty; 
    ligninResults = double.empty;
        
    for i = 1 : r_StdStt     %%% i represents each flux steady-state vector from the static model results %%%
        i       
        p = p_l;

        Xss = steadyStates(i,:)';
        
        StdStt = StdStt_l;
        Vss = StdStt(i,1:38,1);
        Vss(39:42) = StdStt(i,53:56,1);
        Vss(43:56) = StdStt(i,39:52,1);
                
        L_x = squeeze(StdStt(i,1 : metabolite,7:10));
        
        s = 0;
        explored = 0;
        
        solutions = double.empty;
        solutionsAll = double.empty;
        [r_s,~] = size(solutions);
        while ~(s > 0.1 * org_sample && explored > 0.5 * org_sample)% 
            
            [sample,~] = size(p);
            
            for j = 1 : sample     %%% j represents each parameter set from the sample set %%%
                
                results = double.empty;
                explored = org_sample - sample + j;
                if s > 0.1 * org_sample && explored > 0.5 * org_sample % 
                    break
                end
                [i,r_s,explored,s]
                g = zeros(metabolite, n_v + n_e);
                h = zeros(metabolite, n_v + n_e);

                g(1,1) = p(j,1);
                g(3,2) = p(j,2);
                g(5,3) = p(j,3);
                g(6,4) = p(j,4);
                g(9,5) = p(j,5);
                g(11,6) = p(j,6);
                g(12,7) = p(j,7);
                g(13,8) = p(j,8);
                g(14,9) = p(j,9);
                g(15,10) = p(j,10);
                g(12,11) = p(j,11);
                g(13,12) = p(j,12);
                g(17,13) = p(j,13);
                g(18,14) = p(j,14);
                g(18,15) = p(j,15);
                g(19,16) = p(j,16);
                g(20,17) = p(j,17);
                g(21,18) = p(j,18);
                g(21,19) = p(j,19);
                g(23,20) = p(j,20);
                g(24,21) = p(j,21);
                g(25,22) = p(j,22);
                g(26,23) = p(j,23);
                g(27,24) = p(j,24);
                g(2,25) = p(j,25);           
                g(28,26) = p(j,26);
                g(29,27) = p(j,27);
                g(30,28) = p(j,28);
                g(2,29) = p(j,29);
                g(31,30) = p(j,30);
                g(28,31) = p(j,31);
                g(4,32) = p(j,32);
                g(31,33) = p(j,33);
                g(33,34) = p(j,34);
                g(32,35) = p(j,35);
                g(34,36) = p(j,36);
                g(34,37) = p(j,37);
                g(6,38) = p(j,38);
                g(12,39) = p(j,39);
                g(20,40) = p(j,40);
                g(2,41) = p(j,41);
                g(32,42) = p(j,42);
                
                % the parallel pools with less than 100 steady-state
                % concentration
                perc_ER = 5;
                perc_cyt = 25;
                Xss(1) = (100 - perc_ER / (1 - r)) * p(j,43) + perc_ER / (1 - r);% none of the ER pools contain less than perc_ER% of the total pool              
                Xss(17) = (100 - perc_cyt / r) * p(j,44) + perc_cyt / r;% none of the ER pools contain less than perc_cyt% of the total pool
                Xss(24) = (100 - perc_cyt / r) * p(j,45) + perc_cyt / r;
                Xss(25) = (100 - perc_cyt / r) * p(j,46) + perc_cyt / r;
                %%%

                h(10,5) = - p(j,47);
                h(10,6) = - p(j,48);
                h(12,6) = - p(j,49);

                KO = g + h;

                % the parallel pools with higher than 100 steady-state
                % concentration
                Xss(10) = (100 - (1 - r) * Xss(1)) / r;
                Xss(4) = (100 - r * Xss(17)) / (1 - r);
                Xss(7) = (100 - r * Xss(24)) / (1 - r);
                Xss(8) = (100 - r * Xss(25)) / (1 - r);
                %unparallel pools
                Xss(9) = 100 / r;
                Xss(11) = 100 / r;
                Xss(15) = 100;
                Xss(30) = 100;
                Xss(26) = 100 / r;
                Xss(27) = 100 / r;            
                               
                Q1 = [Xss(10) - Xss(1), Xss(2) - Xss(12), Xss(16) - Xss(3), Xss(4) - Xss(17), Xss(22) - Xss(5), ...
                      Xss(23) - Xss(6), Xss(7) - Xss(24), Xss(8) - Xss(25), Xss(28) - Xss(13), Xss(29) - Xss(14), ...
                      Xss(31) - Xss(18), Xss(32) - Xss(20), Xss(33) - Xss(19), Xss(34) - Xss(21)];

                Q2 = prod(repmat(Xss,1,n_v + n_e) .^ KO);

                SS = [Xss .* (ones(metabolite,4) - L_x); Xss .* L_x];

                a_V = Vss(1:42) ./ Q2;           
                a_D = Vss(43:56) ./ Q1;

                %%initial conditions%%
                tf = 1e4;
                tspan = 0 : tf/1000 : tf;
                %%%%% down regulations %%%%%
                CTRL_Phe = 1;
                CTRL_Tyr = 1;
                CA_Phe = 1;
                CA_Tyr = 1;
                PCA500_Phe = 1;
                PCA500_Tyr = 1;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                pass = 1;
                %%
                if CTRL_Phe == 1 && pass == 1
                    pass = 0;
                    bolus = [0, 0];  
                    CTRL = 1;
                    StartTime = clock;
                    options = odeset('Events',@(t, x) events(t, x, pools,CTRL,StartTime));
                    x0 = 1.1 * SS(:,1);
                    A = ([a_V a_D])';
                    [t,x,tE,xE,IE] = ode15s(@(t,x,dx)brach(t,x,A,KO,n_v + n_e,metabolite,I_L(:,1),I_U(:,1),S,S_p,D_xp,D_xn,r,bolus,CA_U, PCA500_U,pools,Xss),tspan,x0,options);           
                    [r_t,~] = size(t);
                    if r_t < 2
                        continue
                    end

                    if tE < tf
                        if IE(length(IE)) ~= 2
                            if ~isempty(solutions)
                                s = s + 1;
                            end
                            continue
                        end
                    end

                    x_L = x(r_t,metabolite + 1 : 2 * metabolite);
                    x_U = x(r_t,1 : metabolite);
                    xf = x(r_t,metabolite + 1 : 2 * metabolite) + x(r_t,1 : metabolite);
                    U = x_U ./ xf;
                    L = x_L ./ xf;
                    
                    if norm(((x(r_t,:))' - SS(:,1))) / norm(SS(:,1)) > 5e-2
                        continue
                    end

                    %%fluxes%%
                    W = prod(repmat(xf',1,n_v + n_e).^ KO)';
                    Vf = A(1:42) .* W;    

                    V_L = S_p * L' .* Vf;  
                    V_U = Vf - V_L;

                    L_total_lignin = r * (V_L(10) + V_L(20) + V_L(24)) + (1 - r) * (V_L(28) + V_L(38));
                    U_total_lignin = r * (V_U(10) + V_U(20) + V_U(24)) + (1 - r) * (V_U(28) + V_U(38));
                    total_lignin = r * (Vf(10) + Vf(20) + Vf(24)) + (1 - r) * (Vf(28) + Vf(38));
                    L_H_lignin = r * (V_L(10)) + (1 - r) * (V_L(28));
                    H_lignin = r * (Vf(10)) + (1 - r) * (Vf(28));
                    L_G_lignin = r * (V_L(20)) + (1 - r) * (V_L(38));
                    G_lignin = r * (Vf(20)) + (1 - r) * (Vf(38));
                    L_S_lignin = r * (V_L(24));
                    S_lignin = r * (Vf(24));
                    SG_ratio = S_lignin / G_lignin;
                    L_pCA = r * V_L(39) + (1 - r) * V_L(41);
                    pCA = r * Vf(39) + (1 - r) * Vf(41);
                    L_FA = r * V_L(40) + (1 - r) * V_L(42);
                    FA = r * Vf(40) + (1 - r) * Vf(42);
                    
                    WT_G_lignin = G_lignin;
                    WT_S_lignin = S_lignin;
                    WT_H_lignin = H_lignin;
                    WT_total_lignin = total_lignin;
                    WT_pCA = pCA;
                    WT_FA = FA;
                    
                    if total_lignin > 0.9 * (I(1) + I(2))
                        if abs(SG_ratio - 1.09) < 0.05
                            if abs((H_lignin / total_lignin - 0.04) / 0.04) < tol
                                if abs((G_lignin / total_lignin - 0.41) / 0.41) < tol
                                    if abs((S_lignin / total_lignin - 0.55) / 0.55) < tol
                                        if abs(L_G_lignin / G_lignin - labeling(1,3)) < 0.05 
                                            if abs(L_S_lignin / S_lignin - labeling(1,4)) < 0.05                                           
                                                if abs(L_total_lignin / total_lignin - labeling(1,5)) < 0.05
                                                    if abs(L_pCA / pCA - labeling(1,6)) < 0.05 && pCA < H_lignin && pCA > 0.4 * H_lignin
                                                        if abs(L_FA / FA - labeling(1,7))  < 0.05 && FA < H_lignin && FA > 0.4 * H_lignin
                                                            if (L_G_lignin / G_lignin > 1.05 * L_S_lignin / S_lignin)
                                                               pass = 1;
                                                               results = [results;L_G_lignin / G_lignin,L_S_lignin / S_lignin,L_total_lignin / total_lignin,L_pCA / pCA,L_FA / FA];
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
                %%
                if CTRL_Tyr == 1 && pass == 1
                    pass = 0;
                    bolus = [0, 0];
                    CTRL = 1;
                    StartTime = clock;
                    options = odeset('Events',@(t, x) events(t, x, pools,CTRL,StartTime));
                    x0 = 1.1 * SS(:,2);
                    A = ([a_V a_D])';
                    [t,x,tE,xE,IE] = ode15s(@(t,x,dx)brach(t,x,A,KO,n_v + n_e,metabolite,I_L(:,2),I_U(:,2),S,S_p,D_xp,D_xn,r,bolus,CA_U, PCA500_U,pools,Xss),tspan,x0,options);           
                    [r_t,~] = size(t);
                    if r_t < 2
                        continue
                    end

                    if tE < tf
                        if IE(length(IE)) ~= 2
                            if ~isempty(solutions)
                                s = s + 1;
                            end
                            continue
                        end
                    end

                    x_L = x(r_t,metabolite + 1 : 2 * metabolite);
                    x_U = x(r_t,1 : metabolite);
                    xf = x_L + x_U;
                    U = x_U ./ xf;
                    L = x_L ./ xf;

                    if norm(((x(r_t,:))' - SS(:,2))) / norm(SS(:,2)) > 5e-2
                        continue
                    end

                    %%fluxes%%
                    W = prod(repmat(xf',1,n_v + n_e).^ KO)';
                    Vf = A(1:42) .* W;    

                    V_L(1:42,1) = S_p * L' .* Vf;  
                    V_U(1:42,1) = Vf - V_L;

                    L_total_lignin = r * (V_L(10) + V_L(20) + V_L(24)) + (1 - r) * (V_L(28) + V_L(38));
                    U_total_lignin = r * (V_U(10) + V_U(20) + V_U(24)) + (1 - r) * (V_U(28) + V_U(38));
                    total_lignin = r * (Vf(10) + Vf(20) + Vf(24)) + (1 - r) * (Vf(28) + Vf(38));
                    L_H_lignin = r * (V_L(10)) + (1 - r) * (V_L(28));
                    H_lignin = r * (Vf(10)) + (1 - r) * (Vf(28));
                    L_G_lignin = r * (V_L(20)) + (1 - r) * (V_L(38));
                    G_lignin = r * (Vf(20)) + (1 - r) * (Vf(38));
                    L_S_lignin = r * (V_L(24));
                    S_lignin = r * (Vf(24));
                    SG_ratio = S_lignin / G_lignin;
                    L_pCA = r * V_L(39) + (1 - r) * V_L(41);
                    pCA = r * Vf(39) + (1 - r) * Vf(41);
                    L_FA = r * V_L(40) + (1 - r) * V_L(42);
                    FA = r * Vf(40) + (1 - r) * Vf(42);
                   
                    if total_lignin > 0.9 * (I(1) + I(2))
                        if abs(SG_ratio - 1.09) < 0.05
                            if abs((H_lignin / total_lignin - 0.04) / 0.04) < tol
                                if abs((G_lignin / total_lignin - 0.41) / 0.41) < tol
                                    if abs((S_lignin / total_lignin - 0.55) / 0.55) < tol
                                        if abs(L_G_lignin / G_lignin - labeling(2,3)) < 0.05 
                                            if abs(L_S_lignin / S_lignin - labeling(2,4)) < 0.05                                           
                                                if abs(L_total_lignin / total_lignin - labeling(2,5)) < 0.05
                                                    if abs(L_pCA / pCA - labeling(2,6)) < 0.05 && pCA < H_lignin && pCA > 0.4 * H_lignin
                                                        if abs(L_FA / FA - labeling(2,7))  < 0.05 && FA < H_lignin && FA > 0.4 * H_lignin
                                                            if (1.05 * L_G_lignin / G_lignin <  L_S_lignin / S_lignin)
                                                               pass = 1;
                                                               results = [results;L_G_lignin / G_lignin,L_S_lignin / S_lignin,L_total_lignin / total_lignin,L_pCA / pCA,L_FA / FA];
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
                %%
                if CA_Phe == 1 && pass == 1
                    pass = 0;
                    bolus = [1, 0];
                    CTRL = 0;
                    StartTime = clock;
                    options = odeset('Events',@(t, x) events(t, x, pools,CTRL,StartTime));
                    x0 = 1.0 * SS(:,1);
                    A = ([a_V a_D])';
                    [t,x,tE,xE,IE] = ode15s(@(t,x,dx)brach(t,x,A,KO,n_v + n_e,metabolite,I_L(:,1),I_U(:,1),S,S_p,D_xp,D_xn,r,bolus,CA_U,PCA500_U,pools,Xss),tspan,x0,options);           
                    [r_t,~] = size(t);
                    if r_t < 2
                        continue
                    end
                    if tE < tf
                        if IE(length(IE)) ~= 2
                            if ~isempty(solutions)
                                s = s + 1;
                            end
                            continue
                        end
                    end

                    x_L = x(r_t,metabolite + 1 : 2 * metabolite);
                    x_U = x(r_t,1 : metabolite);
                    xf = x(r_t,metabolite + 1 : 2 * metabolite) + x(r_t,1 : metabolite);
                    U = x_U ./ xf;
                    L = x_L ./ xf;

                    if norm((x(r_t,:) - x(r_t - 1,:))) / norm(x(r_t,:)) > 1e-2
                        if ~isempty(solutions)
                            s = s + 1;
                        end
                        continue
                    end
                    
                    %%fluxes%%
                    W = prod(repmat(xf',1,n_v + n_e).^ KO)';
                    Vf = A(1:42) .* W(1:42);    

                    V_L = S_p * L' .* Vf;  
                    V_U = Vf - V_L;

                    L_total_lignin = r * (V_L(10) + V_L(20) + V_L(24)) + (1 - r) * (V_L(28) + V_L(38));
                    U_total_lignin = r * (V_U(10) + V_U(20) + V_U(24)) + (1 - r) * (V_U(28) + V_U(38));
                    total_lignin = r * (Vf(10) + Vf(20) + Vf(24)) + (1 - r) * (Vf(28) + Vf(38));
                    L_H_lignin = r * (V_L(10)) + (1 - r) * (V_L(28));
                    H_lignin = r * (Vf(10)) + (1 - r) * (Vf(28));
                    L_G_lignin = r * (V_L(20)) + (1 - r) * (V_L(38));
                    G_lignin = r * (Vf(20)) + (1 - r) * (Vf(38));
                    L_S_lignin = r * (V_L(24));
                    S_lignin = r * (Vf(24));
                    SG_ratio = S_lignin / G_lignin;
                    L_pCA = r * V_L(39) + (1 - r) * V_L(41);
                    pCA = r * Vf(39) + (1 - r) * Vf(41);
                    L_FA = r * V_L(40) + (1 - r) * V_L(42);
                    FA = r * Vf(40) + (1 - r) * Vf(42);

                    if L_G_lignin / G_lignin < 0.01
                        if L_S_lignin / S_lignin < 0.01
                            if L_total_lignin / total_lignin < 0.01
                                if abs(SG_ratio - 1.29) < 0.2
                                   pass = 1;
                                   results = [results;L_G_lignin / G_lignin,L_S_lignin / S_lignin,L_total_lignin / total_lignin,L_pCA / pCA,L_FA / FA];
                                end
                            end
                        end
                    end
                end

                if CA_Tyr == 1 && pass == 1
                    pass = 0;
                    bolus = [1, 0];
                    CTRL = 0;
                    StartTime = clock;
                    options = odeset('Events',@(t, x) events(t, x, pools,CTRL,StartTime));
                    x0 = 1.0 * SS(:,2);
                    A = ([a_V a_D])';
                    [t,x,tE,xE,IE] = ode15s(@(t,x,dx)brach(t,x,A,KO,n_v + n_e,metabolite,I_L(:,2),I_U(:,2),S,S_p,D_xp,D_xn,r,bolus,CA_U,PCA500_U,pools,Xss),tspan,x0,options);           
                    [r_t,~] = size(t);
                    if r_t < 2
                        continue
                    end
                    if tE < tf
                        if IE(length(IE)) ~= 2
                            if ~isempty(solutions)
                                s = s + 1;
                            end
                            continue
                        end
                    end

                    x_L = x(r_t,metabolite + 1 : 2 * metabolite);
                    x_U = x(r_t,1 : metabolite);
                    xf = x(r_t,metabolite + 1 : 2 * metabolite) + x(r_t,1 : metabolite);
                    U = x_U ./ xf;
                    L = x_L ./ xf;

                    if norm((x(r_t,:) - x(r_t - 1,:))) / norm(x(r_t,:)) > 1e-2
                        if ~isempty(solutions)
                            s = s + 1;
                        end
                        continue
                    end

                    %%fluxes%%
                    W = prod(repmat(xf',1,n_v + n_e).^ KO)';
                    Vf = A(1:42) .* W(1:42);    

                    V_L = S_p * L' .* Vf;  
                    V_U = Vf - V_L;

                    L_total_lignin = r * (V_L(10) + V_L(20) + V_L(24)) + (1 - r) * (V_L(28) + V_L(38));
                    U_total_lignin = r * (V_U(10) + V_U(20) + V_U(24)) + (1 - r) * (V_U(28) + V_U(38));
                    total_lignin = r * (Vf(10) + Vf(20) + Vf(24)) + (1 - r) * (Vf(28) + Vf(38));
                    L_H_lignin = r * (V_L(10)) + (1 - r) * (V_L(28));
                    H_lignin = r * (Vf(10)) + (1 - r) * (Vf(28));
                    L_G_lignin = r * (V_L(20)) + (1 - r) * (V_L(38));
                    G_lignin = r * (Vf(20)) + (1 - r) * (Vf(38));
                    L_S_lignin = r * (V_L(24));
                    S_lignin = r * (Vf(24));
                    SG_ratio = S_lignin / G_lignin;
                    L_pCA = r * V_L(39) + (1 - r) * V_L(41);
                    pCA = r * Vf(39) + (1 - r) * Vf(41);
                    L_FA = r * V_L(40) + (1 - r) * V_L(42);
                    FA = r * Vf(40) + (1 - r) * Vf(42);

                    if abs(L_G_lignin / G_lignin - 0.065) < 0.03
                        if abs(L_S_lignin / S_lignin - 0.101) < 0.03
                            if abs(L_total_lignin / total_lignin - 0.081) < 0.03
                                if abs(SG_ratio - 1.29) < 0.2
                                   pass = 1;
                                   results = [results;L_G_lignin / G_lignin,L_S_lignin / S_lignin,L_total_lignin / total_lignin,L_pCA / pCA,L_FA / FA];
                                end
                            end
                        end
                    end             
                end

                if PCA500_Phe == 1 && pass == 1
                    pass = 0;
                    bolus = [0, 1];
                    CTRL = 0;
                    StartTime = clock;
                    options = odeset('Events',@(t, x) events(t, x, pools,CTRL,StartTime));
                    x0 = 1.0 * SS(:,1);
                    A = ([a_V a_D])';
                    [t,x,tE,xE,IE] = ode15s(@(t,x,dx)brach(t,x,A,KO,n_v + n_e,metabolite,I_L(:,1),I_U(:,1),S,S_p,D_xp,D_xn,r,bolus,CA_U,PCA500_U,pools,Xss),tspan,x0,options);           

                    [r_t,~] = size(t);
                    if r_t < 2
                        continue
                    end

                    if tE < tf
                        if IE(length(IE)) ~= 2
                            if ~isempty(solutions)
                                s = s + 1;
                            end
                            continue
                        end
                    end

                    x_L = x(r_t,metabolite + 1 : 2 * metabolite);
                    x_U = x(r_t,1 : metabolite);
                    xf = x(r_t,metabolite + 1 : 2 * metabolite) + x(r_t,1 : metabolite);
                    U = x_U ./ xf;
                    L = x_L ./ xf;

                    if norm((x(r_t,:) - x(r_t - 1,:))) / norm(x(r_t,:)) > 1e-2
                        if ~isempty(solutions)
                            s = s + 1;
                        end
                        continue
                    end

                    %%fluxes%%
                    W = prod(repmat(xf',1,n_v + n_e).^ KO)';
                    Vf = A(1:42) .* W(1:42);    
                    
                    V_L = S_p * L' .* Vf;  
                    V_U = Vf - V_L;

                    L_total_lignin = r * (V_L(10) + V_L(20) + V_L(24)) + (1 - r) * (V_L(28) + V_L(38));
                    U_total_lignin = r * (V_U(10) + V_U(20) + V_U(24)) + (1 - r) * (V_U(28) + V_U(38));
                    total_lignin = r * (Vf(10) + Vf(20) + Vf(24)) + (1 - r) * (Vf(28) + Vf(38));
                    L_H_lignin = r * (V_L(10)) + (1 - r) * (V_L(28));
                    H_lignin = r * (Vf(10)) + (1 - r) * (Vf(28));
                    L_G_lignin = r * (V_L(20)) + (1 - r) * (V_L(38));
                    G_lignin = r * (Vf(20)) + (1 - r) * (Vf(38));
                    L_S_lignin = r * (V_L(24));
                    S_lignin = r * (Vf(24));
                    SG_ratio = S_lignin / G_lignin;
                    L_pCA = r * V_L(39) + (1 - r) * V_L(41);
                    pCA = r * Vf(39) + (1 - r) * Vf(41);
                    L_FA = r * V_L(40) + (1 - r) * V_L(42);
                    FA = r * Vf(40) + (1 - r) * Vf(42);
                    
                    if abs(L_G_lignin / G_lignin - 0.168) < 0.03 
                        if abs(L_S_lignin / S_lignin - 0.158) < 0.03
                            if abs(L_total_lignin / total_lignin - 0.174)< 0.03
                                if abs(SG_ratio - 0.84) < 0.25
                                    pass = 1;
                                    results = [results;L_G_lignin / G_lignin,L_S_lignin / S_lignin,L_total_lignin / total_lignin,L_pCA / pCA,L_FA / FA];
                                    if abs(L_FA / FA - 0.185)  < 0.04
                                        if abs(L_pCA / pCA - 0.125) < 0.04
                                            pass = 1;
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
                if PCA500_Tyr == 1 && pass == 1
                    pass = 0;
                    bolus = [0, 1];
                    CTRL = 0;
                    StartTime = clock;
                    options = odeset('Events',@(t, x) events(t, x, pools,CTRL,StartTime));
                    x0 = 1.0 * SS(:,2);
                    A = ([a_V a_D])';
                    [t,x,tE,xE,IE] = ode15s(@(t,x,dx)brach(t,x,A,KO,n_v + n_e,metabolite,I_L(:,2),I_U(:,2),S,S_p,D_xp,D_xn,r,bolus,CA_U,PCA500_U,pools,Xss),tspan,x0,options);           
                    [r_t,~] = size(t);
                    if r_t < 2
                        continue
                    end

                    if tE < tf
                        if IE(length(IE)) ~= 2
                            if ~isempty(solutions)
                                s = s + 1;
                            end
                            continue
                        end
                    end

                    x_L = x(r_t,metabolite + 1 : 2 * metabolite);
                    x_U = x(r_t,1 : metabolite);
                    xf = x(r_t,metabolite + 1 : 2 * metabolite) + x(r_t,1 : metabolite);
                    U = x_U ./ xf;
                    L = x_L ./ xf;

                    if norm((x(r_t,:) - x(r_t - 1,:))) / norm(x(r_t,:)) > 1e-2
                        if ~isempty(solutions)
                            s = s + 1;
                        end
                        continue
                    end

                    %%fluxes%%
                    W = prod(repmat(xf',1,n_v + n_e).^ KO)';
                    Vf = A(1:42) .* W(1:42);    

                    V_L = S_p * L' .* Vf;  
                    V_U = Vf - V_L;

                    L_total_lignin = r * (V_L(10) + V_L(20) + V_L(24)) + (1 - r) * (V_L(28) + V_L(38));
                    U_total_lignin = r * (V_U(10) + V_U(20) + V_U(24)) + (1 - r) * (V_U(28) + V_U(38));
                    total_lignin = r * (Vf(10) + Vf(20) + Vf(24)) + (1 - r) * (Vf(28) + Vf(38));
                    L_H_lignin = r * (V_L(10)) + (1 - r) * (V_L(28));
                    H_lignin = r * (Vf(10)) + (1 - r) * (Vf(28));
                    L_G_lignin = r * (V_L(20)) + (1 - r) * (V_L(38));
                    G_lignin = r * (Vf(20)) + (1 - r) * (Vf(38));
                    L_S_lignin = r * (V_L(24));
                    S_lignin = r * (Vf(24));
                    SG_ratio = S_lignin / G_lignin;
                    L_pCA = r * V_L(39) + (1 - r) * V_L(41);
                    pCA = r * Vf(39) + (1 - r) * Vf(41);
                    L_FA = r * V_L(40) + (1 - r) * V_L(42);
                    FA = r * Vf(40) + (1 - r) * Vf(42);

                    if abs(L_G_lignin / G_lignin - 0.088) < 0.03
                        if abs(L_S_lignin / S_lignin - 0.11) < 0.03
                            if abs(L_total_lignin / total_lignin - 0.098)< 0.03
                                if abs(SG_ratio - 0.84) < 0.27
                                    pass = 1;
                                    results = [results;L_G_lignin / G_lignin,L_S_lignin / S_lignin,L_total_lignin / total_lignin,L_pCA / pCA,L_FA / FA];
                                    if abs(L_pCA / pCA - 0.065) < 0.05
                                        if abs(L_FA / FA - 0.05)  < 0.04
                                            pass = 1;
                                        end
                                    end
                                end
                            end
                        end
                    end             
                end
                                
                %%
                if pass == 1
                    ligninResults = cat(3,ligninResults,results);
                    solutionsAll = [solutionsAll;p(j,:)];
                    [r_s,~] = size(solutionsAll);
                    V_results = [V_results;StdStt(i,:,:)];
                    X_results = [X_results;Xss'];
                    P_results = [P_results;p(j,:)];
                    if s > 0.1 * org_sample
                        solutions = double.empty;
                    end
                    s = 0;
                    solutions = [solutions;p(j,:)];
                    CG = mean(solutions,1);
                    [r_p,~] = size(p);
                    p_cut = p(j + 1:r_p,:);
                    p = searchSpace(CG,p(j,:),p_cut,org_sample);
                    break
                elseif ~isempty(solutionsAll)
                    s = s + 1;
                end
            end%end of j%
            if j == sample
                break
            end
            
        end%end of while
        if ~isempty(solutionsAll)
            admsbl_pSliced{i} = solutionsAll;
        end
    end%end of i%

    plotligninResults(ligninResults)
    
    for i = 1 : r_StdStt
        solutionsAll = cell2mat(admsbl_pSliced(i));
        [r,c] = size(solutionsAll);
        admsbl_p(1:r,1:c,i) = solutionsAll;
        admsblStdStt(i,1:2) = [i,r]
    end

    if ~isempty(admsbl_p)
        P = solutionIntegrator(admsbl_p);
    end
    
end

function dx = brach(t,x,A,KO,flux,n_x,I_L,I_U,S,S_p,D_xp,D_xn,r,bolus,CA_U,PCA500_U,pools,x_ss)

    rc = r / (1 - r);
    
    dx = zeros(2 * n_x,1);

    L_U = x(1 : n_x) + x(n_x + 1 : 2 * n_x);
    U = x(1 : n_x) ./ L_U;
    L = x(n_x + 1 : 2 * n_x) ./ L_U;   

    W = prod(repmat(L_U,1,flux) .^ KO)';

    V(1:42,1) = A(1:42) .* W(1:42); 

    V_L(1:42,1) = S_p * L .* V;
   
    V_U(1:42,1) = V - V_L;    

    V_U(43:56) = A(43:56) .* (D_xp - D_xn) * x(1 : n_x);
    V_L(43:56) = A(43:56) .* (D_xp - D_xn) * x(n_x + 1 : 2 * n_x);
    
    dx(1:34) = S * V_U + (1 / r) * I_U;
    dx(35:68) = S * V_L + (1 / r) * I_L;

    dx(9) = dx(9) - (L_U(9) > x_ss(9)) * (0.8 * (L_U(9) - x_ss(9)) * U(9));
    dx(43) = dx(43) - (L_U(9) > x_ss(9)) * (0.8 * (L_U(9) - x_ss(9)) * L(9));
    dx(11) = dx(11) - (L_U(11) > x_ss(11)) * (0.8 * (L_U(11) - x_ss(11)) * U(11));
    dx(45) = dx(45) - (L_U(11) > x_ss(11)) * (0.8 * (L_U(11) - x_ss(11)) * L(11));
    
    if bolus(1)
        dx(10) = dx(10) + (1 / r) * CA_U;
    elseif bolus(2)
        dx(12) = dx(12) + (1 / r) * PCA500_U;
    end

end
function [value,isterminal,direction] = events(t,x,pools,CTRL,StartTime)
    L_U = x(1 : 34) + x(35 : 68);
    persistent Y0 T s T0;

    if t == 0
        Y0 = zeros(68,1);
        T0 = t;
        s = 0;
    end
    A = 100;
    if t > T0 + 10
        A = (norm(abs(x - Y0)) / norm(x)) - ((CTRL > 0) * 1e-3 + (CTRL == 0) * 1e-2);
        Y0 = x;
        T0 = t;
    end
    
    if abs(t - T) < 1e-7
        s = s + 1;
    else
        s = 0;
    end
    TimeElapsed = clock - StartTime;
    value = [min(L_U) - 1; A; 3000 - max(pools * L_U); 5000 - s; 10 - TimeElapsed(end)];% pools * L_U gives the total concentration in cell
    if value(1) < 1e-1
        isterminal = [1; 0; 1; 1; 1];
    elseif value(3) < 1e-1
        isterminal = [1; 0; 1; 1; 1];
    elseif value(4) < 1
        isterminal = [1; 0; 1; 1; 1];
    elseif value(5) < 1
        isterminal = [1; 0; 1; 1; 1];
    else
        isterminal = [1; 1; 1; 1; 1];% going to zero or overaccumulation are in priority to abort to convergence (zero derivatives)
    end
    direction = [0; 0; 0; 0; 0];
    T = t;
end
function [] = plotligninResults(ligninResults)
    Labels = {'$$G$$','$$S$$','$$T$$','$$\textit{p}CA$$','$$FA$$'};
    figure('position',[2000 0 1200 1000]);
    ha = tight_subplot(3, 2, 0.05, [0.06 0.04], [0.12 0.03]);
    data = [0.223   0.210  0.222   0.21     0.23   % labelled phe CTRL
            0.165   0.181  0.186   0.17     0.13   % labelled tyr CTRL
            0       0      0       0        0
            0.065   0.101  0.081   0        0
            0.168   0.158  0.174   0.125    0.185
            0.088   0.11   0.098   0.065    0.05];
    err = [3.6	2.9	3.2 1.0 1.0
                  0.5	0.2	1.5 1.0 1.0
                  0     0   0   0   0
                  0.2	0.2	1.0 0   0
                  1.8	2.0	2.0 2.0 1.5
                  1.3	1.2	1.4 1.0 1.0];
    for i = 1 : 6
        axes(ha(i));
        if i == 1 || i == 2 || i == 5 || i == 6
            boxplot(100 * squeeze(ligninResults(i,:,:))','labels',Labels);hold on     
            scatter([1:5],100 * data(i,:),100,'*','MarkerEdgeColor',[0 .5 .5],'LineWidth',1.8);hold on
        elseif i == 3 || i == 4
            boxplot(100 * squeeze(ligninResults(i,1:3,:))','labels',Labels(1:3));hold on     
            scatter([1:3],100 * data(i,1:3),100,'*','MarkerEdgeColor',[0 .5 .5],'LineWidth',1.8);hold on
        end

        h = findobj(gca,'tag','Outliers');
        delete(h)
        h = findobj(gca,'Tag','Box');
        set(h,'linewidth',1.5);
        h = gca;
        set(gca,'fontsize',16,'FontName', 'Times new roman');
        xlim([0.5 5.5]);
        ylim(100 * [-0.01 0.3]);
        grid on;
        if i == 1 title('^1^3C_9- phenylalanine and unlabeled tyrosine'); end
        if i == 2 title('^1^3C_9- tyrosine and unlabeled phenylalanine'); end
        if i == 1 ylabel('Control','FontWeight','bold'); end
        if i == 3 ylabel('100 \muM CA','FontWeight','bold'); end
        if i == 5 ylabel('500 \muM \it p CA','FontWeight','bold'); end
        h.XAxis.TickLabelInterpreter = 'latex';
    end
end

