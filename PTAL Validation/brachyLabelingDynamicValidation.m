function [ ] = brachyLabelingDynamicValidation()
    clc
    
    load BdPTAL_TGSH1_Lphe_Ltyr
            
    StdStt_l = V_val;
    [r_StdStt,~,~] = size(StdStt_l);
    X_l =  X_val;
    p_l =  P_val;
    enzymeProfile = ep;
    [sample,~] = size(enzymeProfile);
      
    s = 0;

    r = 0.90;
    rc = r / (1 - r);
    
    [S,S_p,D_xp,D_xn,pools,S_x] = Stoi(r);
   
    flux = 56;
    metabolite = 34;
    
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
     
    admsbl = zeros(r_StdStt,sample);
    p = p_l;
    StdStt = StdStt_l;

    ligninResults = double.empty;
    
    for j = 1 : r_StdStt
        j
        tic
        
        Xss = X_l(j,:)';        
        
        Vss = StdStt(j,1:38,1);
        Vss(39:42) = StdStt(j,53:56,1);
        Vss(43:56) = StdStt(j,39:52,1);        
        
        L_x = squeeze(StdStt(j,1 : metabolite,7:10));
                        
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
        
        h(10,5) = 2.5 * h(10,5);
        Km5 = 1.3 * Xss(10);
        a_V(5) = Vss(5) * (Km5 ^ h(10,5) + Xss(10) ^ h(10,5)) / (Xss(9) ^ g(9,5) * Xss(10) ^ h(10,5));

        KO = g + h;

        %%initial conditions%%
        tf = 1e4;
        tspan = 0 : tf/1000 : tf;
        %%%%% down regulations %%%%%
        CTRL_Phe = 1;
        BdPTAL_Phe = 1;
        BdPTAL_Tyr = 1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        pass = 1;
        
        if CTRL_Phe == 1 && pass == 1
            pass = 0;
            bolus = [0, 0];  
            CTRL = 1;
            StartTime = clock;
            options = odeset('Events',@(t, x) events(t, x, pools,CTRL,StartTime));
            x0 = 1.1 * SS(:,1);
            A = ([a_V a_D])';

            [t,x,tE,xE,IE] = ode15s(@(t,x,dx)brach(t,x,A,KO,n_v + n_e,metabolite,I_L(:,1),I_U(:,1),S,S_p,D_xp,D_xn,r,bolus,pools,Xss,Km5),tspan,x0,options);           

            [r_t,~] = size(t);
            if r_t < 2
                continue
            end

            if tE < tf
                if IE(length(IE)) ~= 2
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
            Vf(5) = Vf(5) / (Km5 ^ KO(10,5) + xf(10) ^ KO(10,5));
            Vf_C = Vf;

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
        if pass == 1
            
            for i = 1 : sample     %%% i represents each enzyme purturbation vector from the sample set %%%
                 
                results = double.empty;
                if BdPTAL_Phe == 1

                    pass_l = 0;
                    bolus = [0, 0];
                    CTRL = 0;
                    StartTime = clock;
                    options = odeset('Events',@(t, x) events(t, x, pools,CTRL,StartTime));
                    x0 = [Xss; zeros(metabolite,1)];
                    A = ([a_V a_D])';
                    A(5) = A(5) * (1 - enzymeProfile(i,13) * (1 - 0.5));
                    %TAL
                    A(6) = A(6) * (1 - enzymeProfile(i,1) * (1 - 0.5));
                    %C4H
                    A(1) = A(1) * (1 - enzymeProfile(i,2) * (1 - 1.14));
                    %4CL
                    A(7) = A(7) * (1 - (0.4 * enzymeProfile(i,3) + 0.6) * (1 - 0.3));
                    A(14) = A(14) * (1 - (0.4 * enzymeProfile(i,3) + 0.6) * (1 - 0.3));
                    A(17) = A(17) * (1 - (0.4 * enzymeProfile(i,3) + 0.6) * (1 - 0.3));
                    A(25) = A(25) * (1 - (0.4 * enzymeProfile(i,3) + 0.6) * (1 - 0.3));
                    A(33) = A(33) * (1 - (0.4 * enzymeProfile(i,3) + 0.6) * (1 - 0.3));
                    A(35) = A(35) * (1 - (0.4 * enzymeProfile(i,3) + 0.6) * (1 - 0.3));
                    %HCT1
                    A(12) = A(12) * (1 - enzymeProfile(i,4) * (1 - 0.38));
                    A(31) = A(31) * (1 - enzymeProfile(i,4) * (1 - 0.38));
                    %HCT2
                    A(13) = A(13) * (1 - enzymeProfile(i,5) * (1 - 0.87));
                    A(32) = A(32) * (1 - enzymeProfile(i,5) * (1 - 0.87));
                    %C3'H
                    A(2) = A(2) * (1 - enzymeProfile(i,6) * (1 - 1.19));
                    %COMT
                    A(30) = A(30) * (1 - enzymeProfile(i,7) * (1 - 0.89));
                    A(15) = A(15) * (1 - enzymeProfile(i,7) * (1 - 0.89));
                    A(21) = A(21) * (1 - enzymeProfile(i,7) * (1 - 0.89));
                    A(22) = A(22) * (1 - enzymeProfile(i,7) * (1 - 0.89));
                    %CCoAOMT
                    A(16) = A(16) * (1 - enzymeProfile(i,8) * (1 - 0.73));
                    A(34) = A(34) * (1 - enzymeProfile(i,8) * (1 - 0.73));
                    %CCR
                    A(8) = A(8) * (1 - enzymeProfile(i,9) * (1 - 1.33));
                    A(18) = A(18) * (1 - enzymeProfile(i,9) * (1 - 1.33));
                    A(26) = A(26) * (1 - enzymeProfile(i,9) * (1 - 1.33));
                    A(36) = A(36) * (1 - enzymeProfile(i,9) * (1 - 1.33));
                    %CAD
                    A(9) = A(9) * (1 - enzymeProfile(i,10) * (1 - 1.2));
                    A(19) = A(19) * min((1 - enzymeProfile(i,10) * (1 - 1.2)),(1 - enzymeProfile(i,9) * (1 - 1.33)));
                    A(23) = A(23) * (1 - enzymeProfile(i,10) * (1 - 1.2));
                    A(27) = A(27) * (1 - enzymeProfile(i,10) * (1 - 1.2));
                    A(37) = A(37) * min((1 - enzymeProfile(i,10) * (1 - 1.2)),(1 - enzymeProfile(i,9) * (1 - 1.33)));
                    %F5H
                    A(3) = A(3) * (1 - enzymeProfile(i,11) * (1 - 1.19));
                    A(4) = A(4) * (1 - enzymeProfile(i,11) * (1 - 1.19));
                    %C3H
                    A(11) = A(11) * (1 - enzymeProfile(i,12) * (1 - 0.5));
                    A(29) = A(29) * (1 - enzymeProfile(i,12) * (1 - 0.5));
       
                    %decrease in carbon influx
                    I_UU = I_U * (1 - 0.7 * enzymeProfile(i,14) * (1 - 0.54));
                    I_LL = I_L * (1 - 0.7 * enzymeProfile(i,14) * (1 - 0.54));
                    
                    [t,x,tE,xE,IE] = ode15s(@(t,x,dx)brach(t,x,A,KO,n_v + n_e,metabolite,I_LL(:,1),I_UU(:,1),S,S_p,D_xp,D_xn,r,bolus,pools,Xss,Km5),tspan,x0,options);           

                    [r_t,~] = size(t);
                    if r_t < 2
                        continue
                    end

                    if tE < tf
                        if IE(length(IE)) ~= 2
                            continue
                        end
                    end

                    x_L = x(r_t,metabolite + 1 : 2 * metabolite);
                    x_U = x(r_t,1 : metabolite);
                    xf = x(r_t,metabolite + 1 : 2 * metabolite) + x(r_t,1 : metabolite);
                    U = x_U ./ xf;
                    L = x_L ./ xf;

                    if norm((x(r_t,:) - x(r_t - 1,:))) / norm(x(r_t,:)) > 1e-2
                        continue
                    end

                    %%fluxes%%
                    W = prod(repmat(xf',1,n_v + n_e).^ KO)';
                    Vf = A(1:42) .* W(1:42);    
                    Vf(5) = Vf(5) / (Km5 ^ KO(10,5) + xf(10) ^ KO(10,5));

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
                    
                    if abs(total_lignin / WT_total_lignin - 0.54) < 0.1
                        if abs(G_lignin / WT_G_lignin - 0.58) < 0.1
                            if abs(S_lignin / WT_S_lignin - 0.45) < 0.1
                                if abs(H_lignin / WT_H_lignin - 0.84) < 0.1
                                    if abs(L_total_lignin / total_lignin - 0.26)< 0.05
                                        if abs(L_G_lignin / G_lignin - 0.28) < 0.06 
                                            if abs(L_S_lignin / S_lignin - 0.24) < 0.05
                                                if abs(L_H_lignin / H_lignin - 0.30) < 0.06
                                                    pass_l = 1;
                                                    results = [results;G_lignin / WT_G_lignin,S_lignin / WT_S_lignin,H_lignin / WT_H_lignin,total_lignin / WT_total_lignin,L_G_lignin / G_lignin,L_S_lignin / S_lignin,L_H_lignin / H_lignin,L_total_lignin / total_lignin];
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end             
                end
                if BdPTAL_Tyr == 1 && pass_l == 1
                    pass_l = 0;
                    bolus = [0, 0];
                    CTRL = 0;
                    StartTime = clock;
                    options = odeset('Events',@(t, x) events(t, x, pools,CTRL,StartTime));
                    x0 = [Xss; zeros(metabolite,1)];
                    A = ([a_V a_D])';
                    A(5) = A(5) * (1 - enzymeProfile(i,13) * (1 - 0.5));
                    %TAL
                    A(6) = A(6) * (1 - enzymeProfile(i,1) * (1 - 0.5));
                    %C4H
                    A(1) = A(1) * (1 - enzymeProfile(i,2) * (1 - 1.14));
                    %4CL
                    A(7) = A(7) * (1 - (0.4 * enzymeProfile(i,3) + 0.6) * (1 - 0.3));
                    A(14) = A(14) * (1 - (0.4 * enzymeProfile(i,3) + 0.6) * (1 - 0.3));
                    A(17) = A(17) * (1 - (0.4 * enzymeProfile(i,3) + 0.6) * (1 - 0.3));
                    A(25) = A(25) * (1 - (0.4 * enzymeProfile(i,3) + 0.6) * (1 - 0.3));
                    A(33) = A(33) * (1 - (0.4 * enzymeProfile(i,3) + 0.6) * (1 - 0.3));
                    A(35) = A(35) * (1 - (0.4 * enzymeProfile(i,3) + 0.6) * (1 - 0.3));
                    %HCT1
                    A(12) = A(12) * (1 - enzymeProfile(i,4) * (1 - 0.38));
                    A(31) = A(31) * (1 - enzymeProfile(i,4) * (1 - 0.38));
                    %HCT2
                    A(13) = A(13) * (1 - enzymeProfile(i,5) * (1 - 0.87));
                    A(32) = A(32) * (1 - enzymeProfile(i,5) * (1 - 0.87));
                    %C3'H
                    A(2) = A(2) * (1 - enzymeProfile(i,6) * (1 - 1.19));
                    %COMT
                    A(30) = A(30) * (1 - enzymeProfile(i,7) * (1 - 0.89));
                    A(15) = A(15) * (1 - enzymeProfile(i,7) * (1 - 0.89));
                    A(21) = A(21) * (1 - enzymeProfile(i,7) * (1 - 0.89));
                    A(22) = A(22) * (1 - enzymeProfile(i,7) * (1 - 0.89));
                    %CCoAOMT
                    A(16) = A(16) * (1 - enzymeProfile(i,8) * (1 - 0.73));
                    A(34) = A(34) * (1 - enzymeProfile(i,8) * (1 - 0.73));
                    %CCR
                    A(8) = A(8) * (1 - enzymeProfile(i,9) * (1 - 1.33));
                    A(18) = A(18) * (1 - enzymeProfile(i,9) * (1 - 1.33));
                    A(26) = A(26) * (1 - enzymeProfile(i,9) * (1 - 1.33));
                    A(36) = A(36) * (1 - enzymeProfile(i,9) * (1 - 1.33));
                    %CAD
                    A(9) = A(9) * (1 - enzymeProfile(i,10) * (1 - 1.2));
                    A(19) = A(19) * min((1 - enzymeProfile(i,10) * (1 - 1.2)),(1 - enzymeProfile(i,9) * (1 - 1.33)));
                    A(23) = A(23) * (1 - enzymeProfile(i,10) * (1 - 1.2));
                    A(27) = A(27) * (1 - enzymeProfile(i,10) * (1 - 1.2));
                    A(37) = A(37) * min((1 - enzymeProfile(i,10) * (1 - 1.2)),(1 - enzymeProfile(i,9) * (1 - 1.33)));
                    %F5H
                    A(3) = A(3) * (1 - enzymeProfile(i,11) * (1 - 1.19));
                    A(4) = A(4) * (1 - enzymeProfile(i,11) * (1 - 1.19));
                    %C3H
                    A(11) = A(11) * (1 - enzymeProfile(i,12) * (1 - 0.5));
                    A(29) = A(29) * (1 - enzymeProfile(i,12) * (1 - 0.5));
                    %decrease in carbon influx
                    I_UU = I_U * (1 - 0.7 * enzymeProfile(i,14) * (1 - 0.54));
                    I_LL = I_L * (1 - 0.7 * enzymeProfile(i,14) * (1 - 0.54));
                    
                    [t,x,tE,xE,IE] = ode15s(@(t,x,dx)brach(t,x,A,KO,n_v + n_e,metabolite,I_LL(:,2),I_UU(:,2),S,S_p,D_xp,D_xn,r,bolus,pools,Xss,Km5),tspan,x0,options);           
    
                    [r_t,~] = size(t);
                    if r_t < 2
                        continue
                    end

                    if tE < tf
                        if IE(length(IE)) ~= 2
                            continue
                        end
                    end

                    x_L = x(r_t,metabolite + 1 : 2 * metabolite);
                    x_U = x(r_t,1 : metabolite);
                    xf = x(r_t,metabolite + 1 : 2 * metabolite) + x(r_t,1 : metabolite);
                    U = x_U ./ xf;
                    L = x_L ./ xf;

                    if norm((x(r_t,:) - x(r_t - 1,:))) / norm(x(r_t,:)) > 1e-2
                        continue
                    end

                    %%fluxes%%
                    W = prod(repmat(xf',1,n_v + n_e).^ KO)';
                    Vf = A(1:42) .* W(1:42); 
                    Vf(5) = Vf(5) / (Km5 ^ KO(10,5) + xf(10) ^ KO(10,5));

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
                  
                    if abs(total_lignin / WT_total_lignin - 0.54) < 0.1
                       if abs(G_lignin / WT_G_lignin - 0.58) < 0.1
                            if abs(S_lignin / WT_S_lignin - 0.45) < 0.1
                                if abs(H_lignin / WT_H_lignin - 0.84) < 0.1
                                    if abs(L_total_lignin / total_lignin - 0.16)< 0.04
                                        if abs(L_G_lignin / G_lignin - 0.15) < 0.04
                                            if abs(L_S_lignin / S_lignin - 0.17) < 0.05
                                                if abs(L_H_lignin / H_lignin - 0.14) < 0.05
                                                    pass_l = 1;
                                                    results = [results;G_lignin / WT_G_lignin,S_lignin / WT_S_lignin,H_lignin / WT_H_lignin,total_lignin / WT_total_lignin,L_G_lignin / G_lignin,L_S_lignin / S_lignin,L_H_lignin / H_lignin,L_total_lignin / total_lignin];
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end             
                end

                if pass_l == 1
                    admsbl(j,i) = 1;
                    ligninResults = cat(3,ligninResults,results);
                    [j,i]
                end

            end%end of i%
        end
        [j,sum(admsbl(j,:))]
        toc
    end%end of j%
    
%     plotligninResultsValidation(ligninResults)

    [V_val,X_val,P_val,ep,admsbl_val] = filterAdmsbl(admsbl,V_val,X_val,P_val,enzymeProfile);

%     save BdPTAL_TGSH1_Lphe_Ltyr V_val X_val P_val ep admsbl_val
    
end

function dx = brach(t,x,A,KO,flux,n_x,I_L,I_U,S,S_p,D_xp,D_xn,r,bolus,pools,x_ss,Km5)

    rc = r / (1 - r);
    
    dx = zeros(2 * n_x,1);

    L_U = x(1 : n_x) + x(n_x + 1 : 2 * n_x);
    U = x(1 : n_x) ./ L_U;
    L = x(n_x + 1 : 2 * n_x) ./ L_U;   

    W = prod(repmat(L_U,1,flux) .^ KO)';

    V(1:42,1) = A(1:42) .* W(1:42); 
    V(5) = V(5) / (Km5 ^ KO(10,5) + L_U(10) ^ KO(10,5));

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
function [V_val,X_val,P_val,eP,admsbl_val] = filterAdmsbl(admsbl,V_results,X_results,P_results,enzymeProfile)
    eP = enzymeProfile(find(sum(admsbl,1)),:);
    V_val = V_results(find(sum(admsbl,2)),:,:);
    X_val = X_results(find(sum(admsbl,2)),:);
    P_val = P_results(find(sum(admsbl,2)),:);
    admsbl_val = admsbl(:,find(sum(admsbl,1)));
    admsbl_val = admsbl_val(find(sum(admsbl_val,2)),:);
end
function [] = plotligninResultsValidation(ligninResults)
    Labels1 = {'$$G$$','$$S$$','$$H$$','$$T$$'};
    Labels2 = {'$$G$$','$$S$$','$$H$$','$$T$$'};
    figure('position',[2000 0 1500 350]);
    ha = tight_subplot(1, 3, 0.05, [0.15 0.1], [0.08 0.03]);
    data = [0.28    0.24    0.30    0.26   % labelled phe
            0.15    0.17    0.14    0.16   % labelled tyr
            0.58     0.45    0.84   0.54]; % lignin composition
            
    for i = 1 : 3
        axes(ha(i));
        if i == 1 || i == 2
            boxplot(100 * squeeze(ligninResults(i,5:8,:))','labels',Labels1);hold on
            scatter([1:4],100 * data(i,:),100,'*','MarkerEdgeColor',[0 .5 .5],'LineWidth',1.8);hold on
            ylim(100 * [-0.01 0.4]);
        elseif i == 3
            boxplot(100 * squeeze(ligninResults(1,1:4,:))','labels',Labels1);hold on
            scatter([1:4],100 * data(i,:),100,'*','MarkerEdgeColor',[0 .5 .5],'LineWidth',1.8);hold on
            ylim(100 * [0.3 1.0]);
        end
        
        h = findobj(gca,'tag','Outliers');
        delete(h)
        h = findobj(gca,'Tag','Box');
        set(h,'linewidth',1.5);
        h = gca;
        set(gca,'fontsize',16,'FontName', 'Times new roman');
        xlim([0.5 4.5]);
        grid on;
        if i == 1 
            title('^1^3C_9- phenylalanine and unlabeled tyrosine');
            ylabel('% of ^1^3C_9 Incorporated','FontWeight','bold');
        end
        if i == 2 title('^1^3C_9- tyrosine and unlabeled phenylalanine'); end
        if i == 3 ylabel('BdPTAL/Wild Type','FontWeight','bold'); end
        h.XAxis.TickLabelInterpreter = 'latex';
    end
end

