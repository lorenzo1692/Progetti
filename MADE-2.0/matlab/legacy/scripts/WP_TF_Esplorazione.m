clearvars; close all; clc
Increm = 1.000; % Per passare a dimensioni a Tamb (1.003)
iter = 0;
maxdim = 30;
%% Materials
E_jckt = 205;
E_ins = 12;
E_cbl_HTS = 120;
E_cbl_LTS = 0;
E_case = 205;
%% Design points
for bbb = 1:4
    counter = 0;
    B0 = ones(1,10)*(5+bbb); % ones(1,10)*13.25;%        
    R0 = [11.25 8.35 8.0 8.3 8.6];  %  
    wb = ones(1,10)*1.4; % linspace(0.7,0.85,10); %  
    A = ones(1,10)*6.5;
    n_TF = 12;              % TF number
    %% 
    for dp=1:10     
        RTFi(dp) = R0(dp)-R0(dp)/A(dp)-wb(dp);    % [m]
        ripple = 0.01;
        RTFo(dp) = (R0(dp)+R0(dp)/A(dp))*(1+1/ripple)^(1/n_TF);    % [m]
        R_VV = RTFi(dp)*1.055;  % [m]                             
        %% WP type
        WP_SC_type  = 'full_LTS'; % full_LTS - full_HTS - Hybrid__
        shape_cable = 'Rect';     % RIS_ - Rect
        %% Ammissibili acciaio
        S_amm_JT = (667*1e6);
        S_amm_VT = (667*1e6);
        %% Parametri operativi
        Mu_0 = 4e-7*pi;
        theta_TF = 2*pi/n_TF;
        amp_corr_R0 = 1/(R0(dp)/RTFi(dp))*(1+1/((R0(dp)/RTFi(dp))^n_TF-1)+1/((RTFo(dp)/R0(dp))^n_TF-1));
        NI = (2*pi*R0(dp)*B0(dp)/Mu_0)/n_TF*1e-6;            % Total TF current [MA]
        B_PHI_TF = Mu_0*n_TF*NI*1e6/(2*pi*RTFi(dp));         % Max field on TF
        B_PHI_0 = B_PHI_TF*amp_corr_R0;
        corr_B_WP       = 1.08;
        B_PHI_TF        = B_PHI_TF*corr_B_WP;  % Max field on TF correction
        NI  = NI*1e6;                              % Total TF current [A]
        S_VV = 90*1e6; % VV Yeld limit
        Tau_discharge = B0(dp)*NI*n_TF*(R0(dp)/A(dp))^2/(R_VV*S_VV);
        %% WP dimensioning
        CASE_w = 2*RTFi(dp)*tan(theta_TF/2); % Case width
        lateral_w_min(dp) = CASE_w*0.05; %0.05; % Case thickness wedge side
        lateral_w_max(dp) = CASE_w*0.15; %0.10; % Case thickness wedge side
        dr_plasma_side  = 0.050; % [m] Spessore case fronte plasma
        toroidal_gap    = 0.030; 
        GoundIns        = 0.006; % Ground insulation WP
        WP_w_max        = 2*(RTFi(dp)-dr_plasma_side-GoundIns)*tan(theta_TF/2)-2*lateral_w_min(dp); % WP width
        WP_w_min        = 2*(RTFi(dp)-dr_plasma_side-GoundIns)*tan(theta_TF/2)-2*lateral_w_max(dp); % WP width
        % WP_w            = WP_w_max-2*GoundIns;
        %% Dati cavo superconduttore Nb3Sn
        INS_grades = 0.0005; % Spessore isolante tra i Grades
        d_fili = 0.82;               % wire diameter [mm] 
        cos_theta = 0.97;         % Cos(theta): considers that wires are twisted, then the real effective section is grater 
        VF = 0.8;                % Void Fraction CICC 27%
        r_cc = 0.005/2;             % Cooling channel radius
        N_fili = [1500 1440 1350 1296 1200 1152 1080 972 960 900 864 810 768 720 675 648 540 486 420 360 324 300 216 180 162 144]; % Possible triplets combinations
        %% Dati configurazioni da iterare
        Iop_max = 9.0e4;  % Operative current - max
        Iop_min = 1.0e4;  % Operative current - min
        Ntlmax = ceil(NI/Iop_min); 
        Ntlmin = ceil(NI/Iop_max);
        min_size_CICC = 0.02; % min CICC size
        max_size_CICC = 0.10; % max CICC size
        min_JT = 0.002; 
        n_turns_max = ceil((WP_w_max-GoundIns*2)/min_size_CICC); % Maximum feasible layers
        n_turns_min = ceil((WP_w_min-GoundIns*2)/max_size_CICC); % Minimum feasible layers
        max_n_layers = min(ceil((0.5-GoundIns*2)/min_size_CICC),floor(Ntlmax/n_turns_min)); % Maximum feasible layers
        min_n_layers = max(ceil((0.5-GoundIns*2)/max_size_CICC),ceil(Ntlmin/n_turns_max));  % Minimum feasible layers
        n_spire_fsbl = Ntlmin:Ntlmax;
        % layers_comb = min_n_layers:max_n_layers;
        %% PANCAKE THEN:
        if mod(n_turns_min,2)==1
        n_turns_min = n_turns_min-1;
        end
        turns_comb = n_turns_min:2:n_turns_max;                    % Turns combinations
        layers_comb = min_n_layers:min(max_n_layers,maxdim);
        % % % else
        % % % turns_comb = n_turns_min:2:n_turns_max;                    % Turns combinations
        % % % n_grades = min_n_layers:2:max_n_layers;
        %% Definition of the desing layouts - all passible combination 
        % for i = 1:size(n_grades,2)       
        %     t{i}   = nchoosek(turns_comb,n_grades(i));
        %     t_{i}  = turns_comb'.*ones(size(turns_comb,2),n_grades(i));
        %     t{i}   = vertcat(t{i},t_{i});     
        %     t_{i}  = sort(t{i});
        %     aa{i}   = sum(t{i},2)>min(n_spire) & sum(t{i},2)<max(n_spire);
        %     x{i}   = aa{i}.*t{i}; 
        %     ind = find(sum(x{i},2)==0);  
        %     x{i}(ind,:) = []; 
        %     combT{i} = x{i};
        %     n_spire_tot{i} = sum(combT{i},2);
        % end
        %% Definition of the desing layouts - limited (equal turns for each grades)
        combT = cell(1,size(layers_comb,2));
        for i = 1:size(layers_comb,2)       
            t{i}  = turns_comb'.*ones(size(turns_comb,2),layers_comb(i));  
            aaa{i}   = sum(t{i},2)>min(n_spire_fsbl) & sum(t{i},2)<max(n_spire_fsbl);
            x{i}   = aaa{i}.*t{i}; 
            ind = find(sum(x{i},2)==0);  
            x{i}(ind,:) = []; 
            combT{i} = x{i};
            n_spire_tot{i} = sum(combT{i},2);
        end
        a=0;
        for i =1:sum(size(combT,2))
            a=a+size(combT{i},1);
        end
        a_ = a*size(lateral_w_min(dp):0.01:lateral_w_max(dp),2);
        for lateral_w = lateral_w_min(dp):0.01:lateral_w_max(dp)
            WP_w = CASE_w-lateral_w*2;  
            for i = 1:size(combT,2) 
                for j = 1:size(combT{i},1)
                    iter = iter+1;
                    n_turns = zeros(1,maxdim);
                    n_turns(1:size(combT{i}(j,:),2)) = combT{i}(j,:);
                    n_layers = layers_comb(i);
                    n_layers_ = ones(1,layers_comb(i));                                          
                    %% Turns in each grade
                    n_spire_ = zeros(1,n_layers);
                    n_spire_(1,1) = sum(n_turns(1:size(combT{i}(j,:),2)));
                    for var = 2:n_layers 
                        n_spire_(1,var) = n_spire_(1,var-1)-n_turns(1,var-1);
                    end   
                    %% Preallocations
                    clear type_cable S_rm_grades B_grades Ic_sc N_sc N_Cu THS N_sc0 S_Cu_HTS R_SC SC_w ...
                    R_J S_REBCO SC_h A_SC A_J Ri Re Iop B Cond_h Cond_w EQV_Cable_Area JT r_cable
                    Cond_h = zeros(1,maxdim); 
                    Cond_w = zeros(1,maxdim); 
                    EQV_Cable_Area = zeros(1,maxdim); 
                    JT = zeros(1,maxdim); 
                    r_cable = zeros(1,maxdim);
                    for k = 1:maxdim
                        type_cable(1,k) = {'---'};
                    end 
                    S_rm_grades = zeros(1,n_layers);
                    B_layers = zeros(1,n_layers);
                    Ic_sc = zeros(1,n_layers);
                    N_sc = zeros(1,maxdim);
                    N_Cu = zeros(1,maxdim);
                    THS = zeros(1,maxdim);
                    N_sc0 = zeros(1,n_layers);
                    R_SC = zeros(1,n_layers);
                    S_Cu_HTS = zeros(1,maxdim);
                    S_REBCO = zeros(1,maxdim); 
                    SC_w = zeros(1,n_layers);
                    R_J = zeros(1,n_layers);
                    SC_h = zeros(1,n_layers);
                    A_SC = zeros(1,n_layers);
                    A_J = zeros(1,n_layers);
                    Ri = zeros(1,n_layers);
                    Re = zeros(1,n_layers);        
                    %% Operative current definition
                    Iop = ceil(NI/n_spire_(1,1));
                    if Iop < Iop_min || Iop > Iop_max
                        continue
                    end                
                    %% Definition of the magnetic field peak in each grade (linera beahvior from B(Re)= Bmax to B(Ri) = 0                 
                    Re(1,1) = RTFi(dp)-dr_plasma_side-GoundIns; % WP innerl-leg outer radius (non-insulated)
                    B_TF = B_PHI_TF; % (n_TF*n_spire_(1,1)*Iop*Mu_0)/(2*pi*Re(1,1)); 
                    B_layers = B_TF.*(n_spire_./n_spire_(1,1)); % Ottengo B considerando un andamento lineare nel WP
                    xB = [10 12 14 16.6];
                    yI = [1420 800 580 200];
                    p = polyfit(xB,yI,2);
                    clear Ic_sc_HJc
                    for var = 1:n_layers                                                                                
                        B = B_layers(1,var);
                        Ic_sc(1,var) = Ic_Nb3Sn_WST(B,d_fili); % Critical current Nb3Sn  
                        Ic_sc_HJc(1,var) = polyval(p,B);
                    end     
                    N_sc0 = Iop./Ic_sc; % Needed wires in LTS cable for each grade   
                    %% Induttanza TF modello shell    
                    clear L Tau_discharge
                    L = Mu_0*R0(dp)*(n_TF*n_spire_(1,1))^2*(1-sqrt(1-(R0(dp)-RTFi(dp))/R0(dp)))/n_TF*1.1;
                    V_MAX = 10000/2; % [V]
                    Tau_discharge = (L*Iop/V_MAX); % [s] - scarico in gruppi di tre bobine
                    E = 1/2*L*n_TF*Iop^2*1e-6;
%             if   > Tau_discharge
%                 Tau_discharge = Tau_discharge_;
%             end
%             V_each = L*Iop/Tau_discharge_*1e-3*3; % Tensione max singolo TF
%             k = 0.5*log(RTFo(dp)/RTFi(dp));
%             r0 = sqrt(RTFo(dp)*RTFi(dp));
%             L = 2*(0.5*Mu_0*r0*(n_TF*n_spire_(1,1))^2*k^2*(besseli(0,k)+2*besseli(1,k)+besseli(2,k)))/n_TF
                    %% 
                    clear LTS_Bmin check_type_LTS check_type_HTS
                    if min(WP_SC_type == 'full_HTS') ~= 0    
                        LTS_Bmin = 0;
                    elseif min(WP_SC_type == 'full_LTS') ~= 0
                        LTS_Bmin = 15.5;
                    elseif min(WP_SC_type == 'Hybrid__') ~= 0  
                        LTS_Bmin = 15.5;   
                    end
                    %%
                    check_type_LTS = zeros(1,n_layers);
                    check_type_HTS = zeros(1,n_layers);
                    for var = 1:n_layers
                        if B_layers(1,var) < LTS_Bmin && N_sc0(1,var) <= max(N_fili)
                            check_type_LTS(1,var) = var;
                            type_cable(1,var) = {'LTS'};
                        else 
                            check_type_HTS(1,var) = var;
                            type_cable(1,var) = {'HTS'};
                        end
                    end    
                    %% solo LTS
                    if  max(check_type_HTS) > 0 && min(WP_SC_type == 'full_LTS') ~= 0
                        continue
                    end
                    %%
                    var_LTS = min(check_type_LTS(check_type_LTS>0)):max(check_type_LTS); 
                    var = min(check_type_LTS(check_type_LTS>0));
                    if min(check_type_LTS(check_type_LTS>0)) >0 
                        % LTS                 
                        a_ = N_fili > N_sc0(1,var);
                        z_ = N_fili(a_);
                        N_sc(1,var) = min(z_);
                        N_sc(1,var) = N_sc0(1,var);
                        S_Cu(1,var) = N_sc(1,var)*pi*d_fili^2/8; % Sezione di rame non segregato 
                        Cu_cross_sect_in_SC(1,var) = S_Cu(1,var)*1e-6; % Sezione di rame non segregato in [m2]
                        Nb3Sn_cross_section(1,var) = S_Cu(1,var)*1e-6;  % controllare rapporto Cu non Cu; qui si considera 1
                        J_copper = 110; % Densità di corrente critica stabilizzatore [A/mm2]
                        if Iop/S_Cu(1,var) <= J_copper % Se supero valore critico correggo aggiungendo fili di rame
                            N_tot(1,var) = N_sc(1,var);
                            N_Cu(1,var) = 0;
                        else
                            N_Cu(1,var) = ceil(4/(pi*d_fili^2)*((Iop/J_copper)-S_Cu(1,var)));
                            N_tot(1,var) = N_sc(1,var)+N_Cu(1,var);  % Se ho aggiunto fili di rame ricalcolo N° fili cavo 
                        end   
                        Cu_cross_sect_seg(1,var) = pi*N_Cu(1,var)*d_fili^2/4*1e-6; % Sezione di rame segregato 
                        Bp = B_layers(1,var);                           
                        Cus = Cu_cross_sect_seg(1,var); % Sezione di rame segregato 
                        Cuin = Cu_cross_sect_in_SC(1,var); % Sezione di rame non segregato
                        NB3SN = Nb3Sn_cross_section(1,var);
                        THS(1,var) = Ths_TF_Nb3Sn(Iop,Bp,Cus,Cuin,NB3SN,Tau_discharge,L);    
                        while THS(1,var)>250 && N_Cu(1,var)>0 
                            N_Cu(1,var) = N_Cu(1,var)+1;
                            Bp = B_layers(1,var);                           
                            Cus = pi*N_Cu(1,var)*d_fili^2/4*1e-6;  
                            Cuin = Cu_cross_sect_in_SC(1,var);
                            NB3SN = Nb3Sn_cross_section(1,var);
                            THS(1,var) = Ths_TF_Nb3Sn(Iop,Bp,Cus,Cuin,NB3SN,Tau_discharge,L);                        
                        end
                        while THS(1,var)<250 && N_Cu(1,var)>0 && THS(1,var) > 0
                            N_Cu(1,var) = N_Cu(1,var)-1;
                            Bp = B_layers(1,var);                           
                            Cus = pi*N_Cu(1,var)*d_fili^2/4*1e-6;  
                            Cuin = Cu_cross_sect_in_SC(1,var);
                            NB3SN = Nb3Sn_cross_section(1,var);
                            THS(1,var) = Ths_TF_Nb3Sn(Iop,Bp,Cus,Cuin,NB3SN,Tau_discharge,L);                        
                        end
                        N_tot(1,var) = N_sc(1,var)+N_Cu(1,var);
                        EQV_Cable_Area_0       = (N_tot(1,var)*pi*d_fili^2/4/cos_theta/VF+(pi*r_cc^2))*1e-6; % Area eqv cavo LTS, 3 cc 7 mm
                        D_eqv                  = (EQV_Cable_Area_0*4/pi)^0.5; % Diametro cavo equivalente da area eqv
                        A_w                    = ((D_eqv+0.0004)^2-D_eqv^2)*pi/4; % Considereo l'area del wrapping in acciaio di spessore 0.4 mm
                        EQV_Cable_Area(1,var)  =  pi/4*D_eqv^2+A_w; % Correggo area equivalente del cavo LTS
                        %
                        THS(1,var_LTS)    = THS(1,var);
                        N_Cu(1,var_LTS)  = N_Cu(1,var);
                        N_sc(1,var_LTS)    = N_sc(1,var);
                        EQV_Cable_Area(1,var_LTS)    = EQV_Cable_Area(1,var);
                    end
                    %%
                    var_HTS = min(check_type_HTS(check_type_HTS>0)):max(check_type_HTS); 
                    var = min(check_type_HTS(check_type_HTS>0));
                    if  min(check_type_HTS(check_type_HTS>0)) >0              
                        %% HTS                     
                        type_cable(1,var)        = {'HTS'};
                        % J_REBCO                = (0.9000*(B_layers(1,var))^2-54.100*(B_layers(1,var))+1.1455e3)*1e6; % Densità di corrente critica REBCO
                        J_REBCO                 = Jc_REBCO(B_layers(1,var))*1e6;
                        S_REBCO(1,var)           = Iop/J_REBCO;
                        J_copper                 = 90e6; % Densità di corrente critica stabilizzatore [A/m2]
                        S_Cu_HTS(1,var)          = Iop/J_copper;        
                        Bp = B_layers(1,var);                           
                        Cus = S_Cu_HTS(1,var); % Sezione di rame segregato 
                        Cuin = 0; % Sezione di rame non segregato
                        NB3SN = 0;%S_REBCO(1,var);
                        THS(1,var) = Ths_TF_Nb3Sn(Iop,Bp,Cus,Cuin,NB3SN,Tau_discharge,L);    
                        while THS(1,var)>250
                            S_Cu_HTS(1,var) = S_Cu_HTS(1,var)*1.05;
                            Bp = B_layers(1,var);                           
                            Cus = S_Cu_HTS(1,var);  
                            Cuin = 0;
                            NB3SN = 0;%S_REBCO(1,var);
                            THS(1,var) = Ths_TF_Nb3Sn(Iop,Bp,Cus,Cuin,NB3SN,Tau_discharge,L);                        
                        end                  
                        while THS(1,var)<250 && S_Cu_HTS(1,var)>0
                            S_Cu_HTS(1,var) = S_Cu_HTS(1,var)*0.99;
                            Bp = B_layers(1,var);                           
                            Cus = S_Cu_HTS(1,var);  
                            Cuin = 0;
                            NB3SN = 0;%S_REBCO(1,var);
                            THS(1,var) = Ths_TF_Nb3Sn(Iop,Bp,Cus,Cuin,NB3SN,Tau_discharge,L);                        
                        end                                         
                        EQV_Cable_HTS_0          = S_REBCO(1,var)+S_Cu_HTS(1,var)+(pi*r_cc^2); % Area eqv cavo HTS 
                        D_eqv(1,var)             = (EQV_Cable_HTS_0*4/pi)^0.5; % Diametro cavo equivalente da area eqv
                        A_wr(1,var)              = ((D_eqv(1,var)+0.0005)^2-D_eqv(1,var)^2)*pi/4;% Jacket tondo interno (design cavo HTS ENEA)
                        EQV_Cable_Area(1,var)    = EQV_Cable_HTS_0+A_wr(1,var); % Correggo area equivalente del cavo HTS
                        %
                        THS(1,var_HTS)    = THS(1,var);
                        S_Cu_HTS(1,var_HTS)    = S_Cu_HTS(1,var);
                        S_REBCO(1,var_HTS)    = S_REBCO(1,var);
                        EQV_Cable_Area(1,var_HTS)    = EQV_Cable_Area(1,var);               
                    end   
                    %%  Definisco sezioni cavi in ogni layers 
                    k_bf = 0.5*log(RTFo(dp)/RTFi(dp)); % k bending free
                    T_bf = 0.5*(k_bf*n_TF*(n_spire_(1,1)*Iop)^2*Mu_0/(2*pi)); % Hoop tension along TF longitudinal axis
                    Fr = (n_TF*(n_spire_(1,1)*Iop)^2*Mu_0/2)*(1-(1/sqrt(1-(RTFi(dp)/R0(dp))^2)));                       
                    %
                    WP_w0(1,1) = 2*Re(1,1)*tan(theta_TF/2)-lateral_w*2-GoundIns*2; % massimo ingombro toroidale WP
                    S_z_JT =  T_bf/(WP_w0(1,1)^2)/2;
                    if S_z_JT > S_amm_JT
                        continue
                    end
                    n_turns_add = 0;
                    p_rs = B_TF^2/(2*Mu_0);                                 % Magnetic pressure, thin WP
                    for var = 1:n_layers 
                        % Definisco dimensioni interne del cavo CICC per l'i-esimo grade
                        R_SC(1,var)              = (0.005)*Increm; %JT(1,var);   % Impogno raggio curvatura corner cavo = a spessore Jakcet
                        INS(1,var)               = (0.001)*Increm; % Isolante di spira
                        Cond_w(1,var)            = WP_w0(1,1)/n_turns(1,1);
                        if cell2sym(type_cable(1,var)) == 'HTS'  
                            E_cbl = E_cbl_HTS;
                        else
                            E_cbl = E_cbl_LTS;
                        end	
                        if min(shape_cable == 'RIS_') ~= 0  
                            SC_w(1,var)              = 2*sqrt(EQV_Cable_Area(1,var)/pi);
                            SC_h(1,var)              = SC_w(1,var);   % Ottengo altezza cavo SC
                            R_J(1,var)               = R_SC(1,var);   % Raggio di curvatura corner Jacket
                            Cond_h(1,var)            = Cond_w(1,var); % Ottengo l'altezza del cavo RIS_ 
                            JT(1,var)                = (Cond_w(1,var)-2*INS(1,var)-SC_w(1,var))/2;
                            Ke_cavo_rad(1,var) = 2*E_jckt*JT(1,var)/Cond_h(1,var)+2*INS(1,var)*E_ins/Cond_h(1,var)+...
                            +(1/(E_cbl*SC_w(1,var)/SC_h(1,var))+2/(E_jckt*Cond_w(1,var)/JT(1,var))+2/(E_ins*Cond_w(1,var)/INS(1,var)))^-1;
                            Ke_cavo_tor(1,var) = 2*E_jckt*JT(1,var)/Cond_w(1,var)+2*INS(1,var)*E_ins/Cond_w(1,var)+...
                            +(1/(E_cbl*SC_h(1,var)/SC_w(1,var))+2/(E_jckt*JT(1,var)/Cond_w(1,var))+2/(E_ins*INS(1,var)/Cond_w(1,var)))^-1;
                        elseif min(shape_cable == 'Rect') ~= 0  
                            JT(1,var) =  min_JT;
                            S_rm_grades(1,var) = 1e30;
                            while S_rm_grades(1,var) > S_amm_JT/1.3
                                JT(1,var) = JT(1,var)+5e-4;
                                SC_w(1,var)              = Cond_w(1,var)-JT(1,var)*2-INS(1,var)*2; % Ottengo larghezza cavo SC 
                                SC_h(1,var)              = (EQV_Cable_Area(1,var)+(4-pi)*R_SC(1,var)^2)/SC_w(1,var); % Ottengo altezza cavo SC
                                R_J(1,var)               = R_SC(1,var) + JT(1,var); % Raggio di curvatura corner Jacket
                                Cond_h(1,var)            = SC_h(1,var)+JT(1,var)*2+INS(1,var)*2; % Ottengo l'altezza del cavo Rect                        
                                Ke_cavo_rad(1,var) = 2*E_jckt*JT(1,var)/Cond_h(1,var)+2*INS(1,var)*E_ins/Cond_h(1,var)+...
                                +(1/(E_cbl*SC_w(1,var)/SC_h(1,var))+2/(E_jckt*Cond_w(1,var)/JT(1,var))+2/(E_ins*Cond_w(1,var)/INS(1,var)))^-1;
                                Ke_cavo_tor(1,var) = 2*E_jckt*JT(1,var)/Cond_w(1,var)+2*INS(1,var)*E_ins/Cond_w(1,var)+...
                                +(1/(E_cbl*SC_h(1,var)/SC_w(1,var))+2/(E_jckt*JT(1,var)/Cond_w(1,var))+2/(E_ins*INS(1,var)/Cond_w(1,var)))^-1;
                                K_jckt = 2*JT(1,var)/Cond_h(1,var)*E_jckt;
                                dcr_jckt = K_jckt/Ke_cavo_rad(1,var);
                                r_steel = (Cond_w(1,var)-2*INS(1,var))/(2*JT(1,var));
                                S_rm_grades(1,var) = p_rs*r_steel*dcr_jckt+S_z_JT;      % Radial memebrane stress Jacket innermost layer 
                            end  
                        end
                        A_SC(1,var)              = ((Cond_w(1,var)-2*INS(1,var))*(Cond_h(1,var)-2*INS(1,var)))-((4-pi)*R_J(1,var)^2); % Sezione cavo non isolato
                        A_J(1,var)               = A_SC(1,var)-EQV_Cable_Area(1,var); % Sezione di acciaio                                
                        Ri(1,var)                = Re(1,var)-Cond_h(1,var)*n_layers_(1,var); % Raggio interno grade i-esimo del WP
                        Re(1,var+1)              = Ri(1,var)-INS_grades;   % Raggio esterno grade (i+1)-esimo del WP     
                        WP_w0(1,var)             = Cond_w(1,var)*n_turns(1,var); % massimo ingombro toroidale WP
                        check_w = 2*Ri(1,var)*tan(theta_TF/2);
                        while (check_w-(WP_w0(1,var)+GoundIns*2))/2 <= toroidal_gap
                        n_turns(1,var) = n_turns(1,var)-2;     
                        % Definisco dimensioni interne del cavo CICC per l'i-esimo grade
                        R_SC(1,var)              = (0.005)*Increm; %JT(1,var);   % Impogno raggio curvatura corner cavo = a spessore Jakcet
                        INS(1,var)               = (0.001)*Increm; % Isolante di spira
                        Cond_w(1,var)            = WP_w0(1,1)/n_turns(1,1);  
                        if cell2sym(type_cable(1,var)) == 'HTS'  
                        E_cbl = E_cbl_HTS;
                        else
                        E_cbl = E_cbl_LTS;
                        end	
                        if min(shape_cable == 'RIS_') ~= 0  
                            SC_w(1,var)              = 2*sqrt(EQV_Cable_Area(1,var)/pi);
                            SC_h(1,var)              = SC_w(1,var);   % Ottengo altezza cavo SC
                            R_J(1,var)               = R_SC(1,var);   % Raggio di curvatura corner Jacket
                            Cond_h(1,var)            = Cond_w(1,var); % Ottengo l'altezza del cavo RIS_ 
                            JT(1,var)                = (Cond_w(1,var)-2*INS(1,var)-SC_w(1,var))/2;
                            Ke_cavo_rad(1,var) = 2*E_jckt*JT(1,var)/Cond_h(1,var)+2*INS(1,var)*E_ins/Cond_h(1,var)+...
                            +(1/(E_cbl*SC_w(1,var)/SC_h(1,var))+2/(E_jckt*Cond_w(1,var)/JT(1,var))+2/(E_ins*Cond_w(1,var)/INS(1,var)))^-1;
                            Ke_cavo_tor(1,var) = 2*E_jckt*JT(1,var)/Cond_w(1,var)+2*INS(1,var)*E_ins/Cond_w(1,var)+...
                            +(1/(E_cbl*SC_h(1,var)/SC_w(1,var))+2/(E_jckt*JT(1,var)/Cond_w(1,var))+2/(E_ins*INS(1,var)/Cond_w(1,var)))^-1;
                        elseif min(shape_cable == 'Rect') ~= 0  
                            JT(1,var) =  min_JT;
                            S_rm_grades(1,var) = 1e30;
                            while S_rm_grades(1,var) > S_amm_JT/1.3
                                JT(1,var) = JT(1,var)+5e-4;
                                SC_w(1,var)              = Cond_w(1,var)-JT(1,var)*2-INS(1,var)*2; % Ottengo larghezza cavo SC 
                                SC_h(1,var)              = (EQV_Cable_Area(1,var)+(4-pi)*R_SC(1,var)^2)/SC_w(1,var); % Ottengo altezza cavo SC
                                R_J(1,var)               = R_SC(1,var) + JT(1,var); % Raggio di curvatura corner Jacket
                                Cond_h(1,var)            = SC_h(1,var)+JT(1,var)*2+INS(1,var)*2; % Ottengo l'altezza del cavo Rect                        
                                Ke_cavo_rad(1,var) = 2*E_jckt*JT(1,var)/Cond_h(1,var)+2*INS(1,var)*E_ins/Cond_h(1,var)+...
                                +(1/(E_cbl*SC_w(1,var)/SC_h(1,var))+2/(E_jckt*Cond_w(1,var)/JT(1,var))+2/(E_ins*Cond_w(1,var)/INS(1,var)))^-1;
                                Ke_cavo_tor(1,var) = 2*E_jckt*JT(1,var)/Cond_w(1,var)+2*INS(1,var)*E_ins/Cond_w(1,var)+...
                                +(1/(E_cbl*SC_h(1,var)/SC_w(1,var))+2/(E_jckt*JT(1,var)/Cond_w(1,var))+2/(E_ins*INS(1,var)/Cond_w(1,var)))^-1;
                                K_jckt = 2*JT(1,var)/Cond_h(1,var)*E_jckt;
                                dcr_jckt = K_jckt/Ke_cavo_rad(1,var);
                                r_steel = (Cond_w(1,var)-2*INS(1,var))/(2*JT(1,var));
                                S_rm_grades(1,var) = p_rs*r_steel*dcr_jckt+S_z_JT;      % Radial memebrane stress Jacket innermost layer 
                            end  
                        end
                            A_SC(1,var)              = ((Cond_w(1,var)-2*INS(1,var))*(Cond_h(1,var)-2*INS(1,var)))-((4-pi)*R_J(1,var)^2); % Sezione cavo non isolato
                            A_J(1,var)               = A_SC(1,var)-EQV_Cable_Area(1,var); % Sezione di acciaio
                            Ri(1,var)                = Re(1,var)-Cond_h(1,var); % Raggio interno grade i-esimo del WP
                            Re(1,var+1)              = Ri(1,var)-INS_grades;   % Raggio esterno grade (i+1)-esimo del WP 
                            WP_w0(1,var)             = Cond_w(1,var)*n_turns(1,var); % massimo ingombro toroidale WP
                            check_w = 2*Ri(1,var)*tan(theta_TF/2);     
                            n_turns_add = n_spire_(1)-sum(n_turns(1:size(combT{i}(j,:),2)));        
                        end  
                    end
                    %%
                    if n_turns((size(combT{i}(j,:),2))) <= 0
                        continue
                    end  
                    %
                    if  var == n_layers && n_turns_add >= 1
                        n_layers_add = ceil(n_turns_add/n_turns((size(combT{i}(j,:),2))));
                        n_layers = n_layers+n_layers_add;               
                        
                        if  (mod(n_turns((size(combT{i}(j,:),2))),2) == 0 && mod(n_turns_add/n_layers_add,2) == 0) || (mod(n_turns((size(combT{i}(j,:),2))),2) ~= 0 && mod(n_turns_add/n_layers_add,2) ~= 0)       
                            pluss = round(n_turns_add/n_layers_add);
                        else
                            pluss = round(n_turns_add/n_layers_add)+1;
                        end 
                        
                        for ivar=1:n_layers_add
                            n_turns(1,(size(combT{i}(j,:),2))+ivar) = pluss;
                        end
                        var_= var;
                        for var = var_+1:n_layers
                            Cond_w(1,var)            = WP_w0(1,1)/n_turns(1,1);
                            EQV_Cable_Area(1,var)    = EQV_Cable_Area(1,var-1);
                            type_cable(1,var)        = type_cable(1,var-1);
                            % Definisco dimensioni interne del cavo CICC per l'i-esimo grade
                            R_SC(1,var)              = (0.005)*Increm; %JT(1,var);   % Impogno raggio curvatura corner cavo = a spessore Jakcet
                            INS(1,var)               = (0.001)*Increm; % Isolante di spira
                            Cond_w(1,var)            = WP_w0(1,1)/n_turns(1,1);  
                            if cell2sym(type_cable(1,var)) == 'HTS'  
                                E_cbl = E_cbl_HTS;
                            else
                                E_cbl = E_cbl_LTS;
                            end	
                            if min(shape_cable == 'RIS_') ~= 0  
                                SC_w(1,var)              = 2*sqrt(EQV_Cable_Area(1,var)/pi);
                                SC_h(1,var)              = SC_w(1,var);   % Ottengo altezza cavo SC
                                R_J(1,var)               = R_SC(1,var);   % Raggio di curvatura corner Jacket
                                Cond_h(1,var)            = Cond_w(1,var); % Ottengo l'altezza del cavo RIS_ 
                                JT(1,var)                = (Cond_w(1,var)-2*INS(1,var)-SC_w(1,var))/2;
                                Ke_cavo_rad(1,var) = 2*E_jckt*JT(1,var)/Cond_h(1,var)+2*INS(1,var)*E_ins/Cond_h(1,var)+...
                                +(1/(E_cbl*SC_w(1,var)/SC_h(1,var))+2/(E_jckt*Cond_w(1,var)/JT(1,var))+2/(E_ins*Cond_w(1,var)/INS(1,var)))^-1;
                                Ke_cavo_tor(1,var) = 2*E_jckt*JT(1,var)/Cond_w(1,var)+2*INS(1,var)*E_ins/Cond_w(1,var)+...
                                +(1/(E_cbl*SC_h(1,var)/SC_w(1,var))+2/(E_jckt*JT(1,var)/Cond_w(1,var))+2/(E_ins*INS(1,var)/Cond_w(1,var)))^-1;
                            elseif min(shape_cable == 'Rect') ~= 0  
                                JT(1,var) = min_JT;
                                S_rm_grades(1,var) = 1e30;
                                while S_rm_grades(1,var) > S_amm_JT/1.3
                                    JT(1,var) = JT(1,var)+5e-4;
                                    SC_w(1,var)              = Cond_w(1,var)-JT(1,var)*2-INS(1,var)*2; % Ottengo larghezza cavo SC 
                                    SC_h(1,var)              = (EQV_Cable_Area(1,var)+(4-pi)*R_SC(1,var)^2)/SC_w(1,var); % Ottengo altezza cavo SC
                                    R_J(1,var)               = R_SC(1,var) + JT(1,var); % Raggio di curvatura corner Jacket
                                    Cond_h(1,var)            = SC_h(1,var)+JT(1,var)*2+INS(1,var)*2; % Ottengo l'altezza del cavo Rect                        
                                    Ke_cavo_rad(1,var) = 2*E_jckt*JT(1,var)/Cond_h(1,var)+2*INS(1,var)*E_ins/Cond_h(1,var)+...
                                    +(1/(E_cbl*SC_w(1,var)/SC_h(1,var))+2/(E_jckt*Cond_w(1,var)/JT(1,var))+2/(E_ins*Cond_w(1,var)/INS(1,var)))^-1;
                                    Ke_cavo_tor(1,var) = 2*E_jckt*JT(1,var)/Cond_w(1,var)+2*INS(1,var)*E_ins/Cond_w(1,var)+...
                                    +(1/(E_cbl*SC_h(1,var)/SC_w(1,var))+2/(E_jckt*JT(1,var)/Cond_w(1,var))+2/(E_ins*INS(1,var)/Cond_w(1,var)))^-1;
                                    K_jckt = 2*JT(1,var)/Cond_h(1,var)*E_jckt;
                                    dcr_jckt = K_jckt/Ke_cavo_rad(1,var);
                                    r_steel = (Cond_w(1,var)-2*INS(1,var))/(2*JT(1,var));
                                    S_rm_grades(1,var) = p_rs*r_steel*dcr_jckt+S_z_JT;      % Radial memebrane stress Jacket innermost layer 
                                end  
                            end
                            A_SC(1,var)              = ((Cond_w(1,var)-2*INS(1,var))*(Cond_h(1,var)-2*INS(1,var)))-((4-pi)*R_J(1,var)^2); % Sezione cavo non isolato
                            A_J(1,var)               = A_SC(1,var)-EQV_Cable_Area(1,var); % Sezione di acciaio                                
                            Ri(1,var)                = Re(1,var)-Cond_h(1,var); % Raggio interno grade i-esimo del WP
                            Re(1,var+1)              = Ri(1,var)-INS_grades;    % Raggio esterno grade (i+1)-esimo del WP 
                            WP_w0(1,var)             = Cond_w(1,var)*n_turns(1,var); % massimo ingombro toroidale WP
                        end  
                    end  
                    if n_layers>maxdim || n_layers<0
                        continue
                    end
                    %% Ricalcolo B grades
                    n_spire_ = zeros(1,n_layers);
                    n_spire_(1,1) = sum(n_turns(1:n_layers));
                    for var = 2:n_layers 
                        n_spire_(1,var) = n_spire_(1,var-1)-n_turns(1,var);
                    end 
                    B_layers = B_TF.*(n_spire_/n_spire_(1,1));
                    Iop = ceil(NI/n_spire_(1,1));
                    %% Dati WP                                
                    WP_h = Re(1,1)-Ri(1,n_layers);            % Altezza del WP
                    A_WP = sum(Cond_h(1:n_layers).*WP_w0(1:n_layers));  % Area WP considerando un rettangolo senza corner 
                    A_SC_tot = A_SC(1)*n_spire_(1);           % Sezione totale superconduttore 
                    A_J_tot = A_J(1)*n_spire_(1);             % Sezione totale acciaio nei jacket  
                    Ri_ = RTFi(dp);                      % Quota a raggio esterno Case
                    Rj_ = Ri_-WP_h-dr_plasma_side-GoundIns*2; % Quota a raggio interno WP con Ground
                    check_w(1,1:n_layers) = 2*Ri(1,1:n_layers).*tan(theta_TF/2);          % massimo ingombro toroidale Case
                    if min((check_w-(WP_w0(1,1:n_layers)+GoundIns*2))/2) < toroidal_gap
                        continue
                    end
                    %% Cechck geometrico sulle dimensioni ottenute nel cavo
                    r_cable(1:n_layers) = Cond_w(1:n_layers)./Cond_h(1:n_layers);             % Aspect ratio cavi di ogni grades
                    %%
                    if min(SC_w) <= 0.005 || min(r_cable(1:n_layers))< 0.99 || min(JT(1:n_layers)) < min_JT || min(r_cable(1:n_layers))> 2.5
                        continue
                    end
                    %% Primary radial stress (Pm+Pb)
                    p_rs = B_TF^2/(2*Mu_0);                                 % Magnetic pressure, thin WP         
                    param = Cond_w(1,n_layers)/SC_w(1,n_layers);
                    if param > 1 && param < 2.8
                        xxx = [1.087,1.136,1.190,1.250,1.316,1.389,1.471,1.563,1.667,1.786,1.923,2.083,2.273,2.500,2.778]; 
                    if cell2sym(type_cable(1,n_layers)) == 'LTS'
                        E_cbl = E_cbl_LTS;
                        yyy = [1.01,1.03,1.06,1.10,1.16,1.22,1.29,1.36,1.43,1.50,1.57,1.64,1.71,1.78,1.85];  
                    elseif cell2sym(type_cable(1,n_layers)) == 'HTS'
                        E_cbl = E_cbl_HTS;
                        yyy = [1.01,1.02,1.02,1.04,1.06,1.07,1.11,1.13,1.16,1.18,1.22,1.27,1.29,1.29,1.31]; 
                    end
                        p = polyfit(xxx,yyy,5);
                        scf = polyval(p,param);        
                    else
                        scf = 1.5;
                    end
                    var = n_layers;
                    Ke_cavo_rad(1,var) = 2*E_jckt*JT(1,var)/Cond_h(1,var)+2*INS(1,var)*E_ins/Cond_h(1,var)+...
                    +(1/(E_cbl*SC_w(1,var)/SC_h(1,var))+2/(E_jckt*Cond_w(1,var)/JT(1,var))+2/(E_ins*Cond_w(1,var)/INS(1,var)))^-1;
                    Ke_cavo_tor(1,var) = 2*E_jckt*JT(1,var)/Cond_w(1,var)+2*INS(1,var)*E_ins/Cond_w(1,var)+...
                    +(1/(E_cbl*SC_h(1,var)/SC_w(1,var))+2/(E_jckt*JT(1,var)/Cond_w(1,var))+2/(E_ins*INS(1,var)/Cond_w(1,var)))^-1;
                    K_jckt = 2*JT(1,var)/Cond_h(1,var)*E_jckt;
                    dcr_jckt = K_jckt/Ke_cavo_rad(1,var);
                    r_steel = (Cond_w(n_layers)-2*INS(n_layers))/(2*JT(n_layers));
                    S_rm = p_rs*r_steel*scf*dcr_jckt; % Radial memebrane stress Jacket innermost layer      
                    %% Cechk membrane stress on Jacket 
                    if S_rm<S_amm_JT/1.3 && S_rm>0 
                        clear Rk_ DTF S_T_JT S_T_VT Ke_WP_rad Ke_WP_tor K_ps_rad K_ps_tor dcr_vault_tor dcr_vault_rad...
                        K_vault_rad K_vault_tor K_lat_rad K_lat_tor S_T_JT S_T_VT
                        DTF = 0.05; % Vaul width 
                        S_T_JT = 1e30; S_T_VT = 1e30;  
                        S_c_VT = 1e30;  S_rm_JT = 1e30;
                        while S_T_VT > S_amm_VT || S_c_VT > S_amm_VT/1.3 || S_T_JT > S_amm_JT || S_rm_JT > S_amm_JT/1.3
                            DTF = DTF+0.001; % Se supero Tresca incremento spessore naso TF
                            Rk_ = Rj_-DTF; % Innermost Case radius 
                            CASE_w_l = 2*Rk_*tan(theta_TF/2); % Case low part width 
                            A_CASE = (2*Ri_*tan(theta_TF/2)+2*Rk_*tan(theta_TF/2))*(Ri_-Rk_)/2;
                            A_VT = A_CASE-A_WP; % Vault section
                            % Rigidezze    
                            Ke_WP_rad = sum(1./(Ke_cavo_rad(1:n_layers).*n_turns(1:n_layers)))^-1;    
                            Ke_WP_tor = sum((n_turns(1:n_layers)./Ke_cavo_tor(1:n_layers)).^-1);   
                            
                            K_ps_rad = E_case/dr_plasma_side*(2*Ri_*tan(theta_TF/2));
                            K_ps_tor = (E_case*dr_plasma_side/(2*Ri_*tan(theta_TF/2)));
                            K_lat_rad = 2*E_case*lateral_w/WP_h;
                            K_lat_tor = (2/(E_case*WP_h/lateral_w))^-1;
                            K_vault_rad = (E_case*Rk_*2*tan(theta_TF/2)/DTF);
                            K_vault_tor = (E_case*DTF/(Rk_*2*tan(theta_TF/2)));
                            
                            %%% Ke_case_rad = (1/(Ke_WP_rad+2*K_lat_rad)+1/K_vault_rad)^-1;
                            Ke_case_rad = (1/Ke_WP_rad+1/K_vault_rad)^-1;
                            Ke_case_tor = K_ps_tor+(1/Ke_WP_tor+1/K_lat_tor)^-1+K_vault_tor;
                            
                            dcr_vault_tor = K_vault_tor/Ke_case_tor;
                            dcr_vault_rad = Ke_WP_rad/Ke_case_rad; 
                            
                            if dcr_vault_rad > 1
                                dcr_vault_rad = 1;
                            end
                            
                            S_rm_JT = S_rm*dcr_vault_rad; % Radial memebrane stress Jacket innermost layer correction                       
                            
                            beta = Rk_/Rj_;
                            S_r_VT = 0; % @ Rk_
                            S_c_VT = 2/(1-beta^2)*p_rs*dcr_vault_tor; % @ Rk_
                            S_l_VT = -1/(1-beta^2)*p_rs; % @ Rk_    
                            
                            k_bf = 0.5*log(RTFo(dp)/RTFi(dp)); % k bending free
                            T_bf = 0.5*(k_bf*n_TF*(n_spire_(1,1)*Iop)^2*Mu_0/(2*pi)); % Hoop tension along TF longitudinal axis
                            
                            S_z =  T_bf/(A_J_tot+A_VT)*1.1; % Normal tension on Case stell section 
                            S_T_VT = (S_z+S_c_VT); % Vault Tresca stress 
                            S_T_JT = S_z+S_rm_JT; % Jacket Tresca stress      
                        end   
                        %%
                        if S_T_VT < S_amm_JT && S_T_JT < S_amm_JT && S_T_JT>0 && S_T_VT>0 && Rk_>0.55
                            counter =  counter+1;  
                            n_cond =  n_spire_(1,1);
                            WP_w = WP_w0(1,1);
                            JENG = Iop/(min(Cond_w(1:n_layers))*min(Cond_h(1:n_layers)))*1e-6;
                            E = 0.5*L*Iop^2*1e-6;                   
                            S_T_VT = S_T_VT*1e-6;
                            S_T_JT = S_T_JT*1e-6;
                            Nose = Rj_-Rk_;
                            R_0 = R0(dp);
                            radial_build = Ri_-Rk_;
                            DATA{bbb}(counter,:) = table(dp,S_T_VT,S_T_JT,R_0,B_PHI_0,B_TF,Iop,JENG,L,E,Ri_,Rj_,Rk_,radial_build,Nose,WP_h,WP_w,...
                            lateral_w,n_cond,n_layers,n_turns,type_cable,Cond_w,Cond_h,JT,r_cable,N_sc,N_Cu,EQV_Cable_Area,S_REBCO,S_Cu_HTS,THS,Tau_discharge);                                                                                                  
                        end             
                    end        
                end  
            end
        end
    end
    %% Scrivo matrice soluzione prodotte
    if size(DATA{bbb},1) > 1
        TitleName = sprintf('TF exploration A%G B0=%G T WP %s.xlsx',A(dp)*10,B0(dp),WP_SC_type);
        writetable(DATA{bbb},TitleName) 
    end
    %%
    figure('units','normalized','outerposition',[0.25 0.1 0.5 0.9])
    Bb = table2array((DATA{bbb}(2:end,6)));
    RB = table2array((DATA{bbb}(2:end,14)));
    Rii = table2array((DATA{bbb}(2:end,11)));
    NIi = table2array((DATA{bbb}(2:end,8)));
    scatter(Rii,RB,20,NIi,'filled') 
    [x1,S] = polyfit(Rii,RB,2);
    [y1,delta] =  polyval(x1,min(Rii):1e-3:max(Rii),S);
    hold on
    % plot(min(Rii):1e-3:max(Rii),y1,'k','linewidth',1)
    colormap(jet);
    set(gca,'FontSize',20);
    colorbar
    c = colorbar;
    c.Label.FontSize = 20;
    c.Label.Color = 'k';
    c.Label.Rotation = 90;
    c.Label.String = 'Jeng [$A/mm^2$]';
    c.Label.Interpreter = 'latex';
    c.TickLabelInterpreter = 'latex';
    mod_max = max(NIi);
    mod_min = min(NIi);
    clim([mod_min mod_max]);
    MAXMIN=get(c,'Limits');
    T = linspace(MAXMIN(1),MAXMIN(2),6);
    set(c,'Ticks',T)
    ylabel('Radial Build TF inner-leg [$m$]','fontsize',20,'Interpreter','latex')
    xlabel('Ri TF Inner-leg [$m$]','fontsize',20,'Interpreter','latex')
    TitleName = sprintf('TF A%G B0=%G T WP %s',A(dp)*10,B0(dp),WP_SC_type);
    title(TitleName,'fontsize',20,'Interpreter','latex')
    grid on
    print('-dpng','-r300',TitleName);
    %%
    figure('units','normalized','outerposition',[0.25 0.1 0.5 0.9])
    wb = table2array((DATA{bbb}(2:end,4)))-table2array((DATA{bbb}(2:end,4)))/A(dp)-table2array((DATA{bbb}(2:end,11)));
    scatter(wb,RB,20,Bb,'filled') 
    [x1,S] = polyfit(wb,RB,2);
    [y1] =  polyval(x1,min(wb):1e-3:max(wb),S);
    hold on
    % plot(min(wb):1e-3:max(wb),y1,'k','linewidth',1)
    colormap(jet);
    set(gca,'FontSize',20);
    colorbar
    c = colorbar;
    c.Label.FontSize = 20;
    c.Label.Color = 'k';
    c.Label.Rotation = 90;
    c.Label.String = '$B_{TF}$ [$T$]';
    c.Label.Interpreter = 'latex';
    c.TickLabelInterpreter = 'latex';
    mod_max = max(Bb);
    mod_min = min(Bb);
    clim([mod_min mod_max]);
    MAXMIN=get(c,'Limits');
    T = linspace(MAXMIN(1),MAXMIN(2),6);
    set(c,'Ticks',T)
    ylabel('Radial Build TF inner-leg [$m$]','fontsize',20,'Interpreter','latex')
    xlabel('VV+SB [$m$]','fontsize',20,'Interpreter','latex')
    TitleName = sprintf('TF A%G B0=%G T WP %s  ',A(dp)*10,B0(dp),WP_SC_type);
    title(TitleName,'fontsize',20,'Interpreter','latex')
    grid on
    print('-dpng','-r300',TitleName);
    %%
    figure('units','normalized','outerposition',[0.25 0.1 0.5 0.9])
    Rkk = table2array((DATA{bbb}(2:end,13)));
    wb = table2array((DATA{bbb}(2:end,4)))-table2array((DATA{bbb}(2:end,4)))/A(dp)-table2array((DATA{bbb}(2:end,11)));
    scatter(wb,Rkk,20,Bb,'filled') 
    [x1,S] = polyfit(wb,Rkk,2);
    [y1] =  polyval(x1,min(wb):1e-3:max(wb),S);
    hold on
    % plot(min(wb):1e-3:max(wb),y1,'k','linewidth',1)
    colormap(jet);
    set(gca,'FontSize',20);
    colorbar
    c = colorbar;
    c.Label.FontSize = 20;
    c.Label.Color = 'k';
    c.Label.Rotation = 90;
    c.Label.String = '$B_{TF}$ [$T$]';
    c.Label.Interpreter = 'latex';
    c.TickLabelInterpreter = 'latex';
    mod_max = max(Bb);
    mod_min = min(Bb);
    clim([mod_min mod_max]);
    MAXMIN=get(c,'Limits');
    T = linspace(MAXMIN(1),MAXMIN(2),6);
    set(c,'Ticks',T)
    ylabel('Rk TF inner-leg [$m$]','fontsize',20,'Interpreter','latex')
    xlabel('VV+SB [$m$]','fontsize',20,'Interpreter','latex')
    TitleName = sprintf('TF A%G B0=%G T WP %s ',A(dp)*10,B0(dp),WP_SC_type);
    title(TitleName,'fontsize',20,'Interpreter','latex')
    grid on
    print('-dpng','-r300',TitleName);
end