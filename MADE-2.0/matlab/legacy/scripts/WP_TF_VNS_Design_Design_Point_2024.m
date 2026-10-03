clearvars; close all; clc
%%
Increm = 1.000; % Per passare a dimensioni a Tamb (1.003)
counter = 0;
iter = 0;
maxdim = 30;

%% Materials

E_jckt = 205;
E_cbl_HTS = 120;
E_cbl_LTS = 0.1;
E_case = 205;
E_ins = 12;

%% Design point
dp = 1;   
WP_SC_type  = 100;      % full_LTS = 100 - full_HTS = 101 - Hybrid = 102
shape_cable = 201;      % RIS = 200 - Rect = 201

%% Input esplorazione
n_TF = 12;  
B0(dp) = 5.6;       
R0(dp) = 2.67;  
A(dp) = 4.25;
wb = 0.86;
RTFi(dp) = R0(dp)-R0(dp)/A(dp)-wb;    % [m]
ripple = 0.01; 
RTFo(dp) = (R0(dp)+R0(dp)/A(dp))*(1/ripple)^(1/n_TF);    % [m]
R_VV = RTFi(dp)*1.05;                  % [m]

%% Storage name
TitleName = sprintf('TF VNS A%G - LTS set 2025.xlsx',A(dp)*10);

%% Ammissibili acciaio
S_amm_JT = (667*1e6);
S_amm_VT = (667*1e6);
S_amm_Cu = (100*1e6);

%% Maximum Tension
V_MAX = 2500; % [V]

%% Parametri operativi
Mu_0 = 4e-7*pi;
theta_TF = 2*pi/n_TF;
amp_corr_R0 = 1/(R0(dp)/RTFi(dp))*(1+1/((R0(dp)/RTFi(dp))^n_TF-1)+1/((RTFo(dp)/R0(dp))^n_TF-1));
NI = (2*pi*R0(dp)*B0(dp)/Mu_0)/n_TF*1e-6;            % Total TF current [MA]

dr_plasma_side  = 0.020; % [m] Spessore case fronte plasma

B_PHI_TF = Mu_0*n_TF*NI*1e6/(2*pi*RTFi(dp)-dr_plasma_side);         % Max field on TF
B_PHI_0 = B_PHI_TF*amp_corr_R0;
NI = NI*1e6;                              % Total TF current [A]
R_TF_Outerleg = RTFo(dp);                 % Outer-leg inner radius
R_TF_Innerleg = RTFi(dp);                 % Inner-leg outer radius
S_VV = 90*1e6; % VV Yeld limit
Tau_discharge1 = B0(dp)*NI*n_TF*(R0(dp)/A(dp))^2/(R_VV*S_VV);
corr_B_WP       = 1.06;
B_PHI_TF        = B_PHI_TF*corr_B_WP;  % Max field on TF correction

%% WP dimensioning
CASE_w          = 2*R_TF_Innerleg*tan(theta_TF/2); % Case width
lateral_w_min   = CASE_w*0.05; % Case thickness wedge side
lateral_w_max   = CASE_w*0.3;

toroidal_gap    = 0.015; 
GoundIns        = 0.005; % Ground insulation WP
WP_w_max        = 2*(R_TF_Innerleg-dr_plasma_side-GoundIns)*tan(theta_TF/2)-2*lateral_w_min(dp); % WP width
WP_w_min        = 2*(R_TF_Innerleg-dr_plasma_side-GoundIns)*tan(theta_TF/2)-2*lateral_w_max(dp); % WP width
% WP_w            = WP_w_max-2*GoundIns;
INS_grades      = 0.0005;         % Spessore isolante tra i Grades

%% Dati configurazioni da iterare
n_grades = 2;
Iop_max = 7.0e4;  % Operative current - max
Iop_min = 1.0e4;  % Operative current - min
Ntlmax = ceil(NI/Iop_min); 
Ntlmin = ceil(NI/Iop_max);
min_size_CICC = 0.015; % min CICC size
max_size_CICC = 0.06; % max CICC size
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
% % % turns_comb = n_turns_min:2:n_turns_max;              % Turns combinations
% % % n_grades = min_n_layers:2:max_n_layers;

%% Definition of the desing layouts - all possible combination
% for i = 1:size(n_grades,2) t{i} = nchoosek(turns_comb,n_grades(i)); t_{i} 
% = turns_comb'.*ones(size(turns_comb,2),n_grades(i)); t{i} = vertcat(t{i},t_{i}); 
% t_{i} = sort(t{i}); aa{i} = sum(t{i},2)>min(n_spire) & sum(t{i},2)<max(n_spire); 
% x{i} = aa{i}.*t{i}; ind = find(sum(x{i},2)==0); x{i}(ind,:) = []; combT{i} = 
% x{i}; n_spire_tot{i} = sum(combT{i},2); end

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
a_ = a*size(lateral_w_min(dp):0.02:lateral_w_max(dp),2);

%% Play each layers-turns combination
for lateral_w = lateral_w_min(dp):0.02:lateral_w_max(dp)
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
            clear type_cable S_rm_grades B_grades Ic_sc N_Sc N_Cu THS N_sc0 S_Cu_HTS r_SC SC_w ...
            R_J S_REBCO SC_h S_CICC S_JT Ri Re Iop B Cond_h Cond_w S_Cable JT r_cable
            Cond_h = zeros(1,maxdim);             Cond_w = zeros(1,maxdim); 
            S_Cable = zeros(1,maxdim);            JT = zeros(1,maxdim); 
            r_cable = zeros(1,maxdim);            tins = zeros(1,maxdim);           
            type_cable = cell(1, maxdim);
            for k = 1:maxdim
                type_cable(1,k) = {'---'};
            end 
            S_rm_grades = zeros(1,n_layers);
            B_layers = zeros(1,n_layers);         Ic_sc = zeros(1,n_layers);
            N_Sc = zeros(1,maxdim);               N_Cu = zeros(1,maxdim);
            THS = zeros(1,maxdim);                N_sc0 = zeros(1,n_layers);
            r_SC = zeros(1,n_layers);             N_tot = zeros(1,n_layers);
            S_Cu_HTS = zeros(1,maxdim);           S_REBCO = zeros(1,maxdim); 
            SC_w = zeros(1,n_layers);             R_J = zeros(1,n_layers);
            SC_h = zeros(1,n_layers);
            S_CICC = zeros(1,n_layers);           S_JT = zeros(1,n_layers);
            Ri = zeros(1,n_layers);               Re = zeros(1,n_layers);

%% Operative current definition
            Iop = ceil(NI/n_spire_(1,1));
            if Iop < Iop_min || Iop > Iop_max
              continue
            end

%% Definition of the magnetic field peak in each grade (linera beahvior from B(Re)= Bmax to B(Ri) = 0
            Re(1,1) = R_TF_Innerleg-dr_plasma_side-GoundIns; % WP innerl-leg outer radius (non-insulated)
            B_TF = B_PHI_TF; % (n_TF*n_spire_(1,1)*Iop*Mu_0)/(2*pi*Re(1,1)); 
            B_layers = B_TF.*(n_spire_./n_spire_(1,1)); % Ottengo B considerando un andamento lineare nel WP          

%% Induttanza TF modello shell
            clear L Tau_discharge
            % L = Mu_0*R0(dp)*(n_TF*n_spire_(1,1))^2*(1-sqrt(1-(R0(dp)-RTFi(dp))/R0(dp)))/n_TF*1.1;
            % L = Mu_0*R0(dp)*(n_TF*n_spire_(1,1))^2*(1-sqrt(1-(RTFo(dp)-RTFi(dp))/2/R0(dp)))/n_TF;

            k_bf = 0.5*log(RTFo(dp)/RTFi(dp));  
            r_bf = sqrt(RTFo(dp)*RTFi(dp));
            Ntot = (n_TF*n_spire_(1,1));
            L = Mu_0*(Ntot*k_bf)^2*r_bf/2*(besseli(0, k_bf)+ 2*besseli(1, k_bf)+ besseli(2, k_bf))/n_TF;

            Tau_discharge2 = (L*Iop/V_MAX); % [s] - scarico in gruppi di n bobine
            Tau_discharge = max([Tau_discharge1 Tau_discharge2 4]);
            E = 1/2*L*n_TF*Iop^2*1e-6;
%%
            for var = 1:n_grades
                B = B_layers(1,var);       
                [type_cable(1,var),N_Cu(1,var),N_Sc(1,var),N_tot(1,var),S_Cable(1,var),S_REBCO(1,var),S_Cu_HTS(1,var),THS(1,var)] ...
                = CICC(B,Iop,Tau_discharge,WP_SC_type); 
            
                type_cable(1,var:n_layers) = type_cable(1,var);
                N_Cu(1,var:n_layers) = N_Cu(1,var);
                N_Sc(1,var:n_layers) = N_Sc(1,var);
                N_tot(1,var:n_layers) = N_tot(1,var);
                S_Cable(1,var:n_layers) = S_Cable(1,var);
                S_REBCO(1,var:n_layers) = S_REBCO(1,var);
                S_Cu_HTS(1,var:n_layers) = S_Cu_HTS(1,var);
                THS(1,var:n_layers) = THS(1,var);
            end

%% Definisco sezioni cavi in ogni layers

            R_WP_IL = R_TF_Innerleg;
            R_WP_OL = R_TF_Outerleg;
            k_bf = 0.5*log(R_WP_OL/R_WP_IL); % k bending free
            T_bf = 0.5*(k_bf*n_TF*(n_spire_(1,1)*Iop)^2*Mu_0/(2*pi)); % Hoop tension along TF longitudinal axis
            Fr = (n_TF*(n_spire_(1,1)*Iop)^2*Mu_0/2)*(1-(1/sqrt(1-(R_WP_IL/R0(dp))^2)));                       
            %
            WP_w0(1,1) = 2*Re(1,1)*tan(theta_TF/2)-lateral_w*2-GoundIns*2; % massimo ingombro toroidale WP
            S_z_JT =  T_bf/(WP_w0(1,1)^2)/2;
            n_turns_add = 0;
            p_rs = B_TF^2/(2*Mu_0);                                 % Magnetic pressure, thin WP
            for var = 1:n_layers 
                % Definisco dimensioni interne del cavo CICC per l'i-esimo grade
                r_SC(1,var)              = (0.005)*Increm; %JT(1,var);   % Impogno raggio curvatura corner cavo = a spessore Jakcet
                tins(1,var)               = (0.001)*Increm; % Isolante di spira
                Cond_w(1,var)            = WP_w0(1,1)/n_turns(1,1);
                if cell2sym(type_cable(1,var)) == 'HTS'  
                    E_cbl = E_cbl_HTS;
                else
                    E_cbl = E_cbl_LTS;
                end	
                if shape_cable == 200
                    SC_w(1,var)              = 2*sqrt(S_Cable(1,var)/pi);
                    SC_h(1,var)              = SC_w(1,var);   % Ottengo altezza cavo SC
                    R_J(1,var)               = r_SC(1,var);   % Raggio di curvatura corner Jacket
                    Cond_h(1,var)            = Cond_w(1,var); % Ottengo l'altezza del cavo RIS_ 
                    JT(1,var)                = (Cond_w(1,var)-2*tins(1,var)-SC_w(1,var))/2;
                    Ke_cavo_rad(1,var)       = 2*E_jckt*JT(1,var)/Cond_h(1,var)+2*tins(1,var)*E_ins/Cond_h(1,var)+...
                    +(1/(E_cbl*SC_w(1,var)/SC_h(1,var))+2/(E_jckt*Cond_w(1,var)/JT(1,var))+2/(E_ins*Cond_w(1,var)/tins(1,var)))^-1;
                    Ke_cavo_tor(1,var)       = 2*E_jckt*JT(1,var)/Cond_w(1,var)+2*tins(1,var)*E_ins/Cond_w(1,var)+...
                    +(1/(E_cbl*SC_h(1,var)/SC_w(1,var))+2/(E_jckt*JT(1,var)/Cond_w(1,var))+2/(E_ins*tins(1,var)/Cond_w(1,var)))^-1;
                elseif shape_cable == 201
                    JT(1,var) =  min_JT;
                    S_rm_grades(1,var) = 1e30;
                    while S_rm_grades(1,var) > S_amm_JT/1.3
                        JT(1,var) = JT(1,var) + 1e-4;
                        SC_w(1,var)              = Cond_w(1,var)-JT(1,var)*2-tins(1,var)*2; % Ottengo larghezza cavo SC 
                        SC_h(1,var)              = (S_Cable(1,var)+(4-pi)*r_SC(1,var)^2)/SC_w(1,var); % Ottengo altezza cavo SC
                        R_J(1,var)               = r_SC(1,var) + JT(1,var); % Raggio di curvatura corner Jacket
                        Cond_h(1,var)            = SC_h(1,var)+JT(1,var)*2+tins(1,var)*2; % Ottengo l'altezza del cavo Rect                        
                        Ke_cavo_rad(1,var)       = 2*E_jckt*JT(1,var)/Cond_h(1,var)+2*tins(1,var)*E_ins/Cond_h(1,var)+...
                        +(1/(E_cbl*SC_w(1,var)/SC_h(1,var))+2/(E_jckt*Cond_w(1,var)/JT(1,var))+2/(E_ins*Cond_w(1,var)/tins(1,var)))^-1;
                        Ke_cavo_tor(1,var)       = 2*E_jckt*JT(1,var)/Cond_w(1,var)+2*tins(1,var)*E_ins/Cond_w(1,var)+...
                        +(1/(E_cbl*SC_h(1,var)/SC_w(1,var))+2/(E_jckt*JT(1,var)/Cond_w(1,var))+2/(E_ins*tins(1,var)/Cond_w(1,var)))^-1;
                        K_jckt                   = 2*JT(1,var)/Cond_h(1,var)*E_jckt;
                        dcr_jckt                 = K_jckt/Ke_cavo_rad(1,var);
                        r_steel                  = (Cond_w(1,var)-2*tins(1,var))/(2*JT(1,var));
                        S_rm_grades(1,var)       = p_rs*r_steel*dcr_jckt+S_z_JT;      % Radial memebrane stress Jacket innermost layer 
                    end  
                end
                S_CICC(1,var)             = ((Cond_w(1,var)-2*tins(1,var))*(Cond_h(1,var)-2*tins(1,var)))-((4-pi)*R_J(1,var)^2); % Sezione cavo non isolato
                S_JT(1,var)               = S_CICC(1,var)-S_Cable(1,var); % Sezione di acciaio                                
                Ri(1,var)                = Re(1,var)-Cond_h(1,var)*n_layers_(1,var); % Raggio interno grade i-esimo del WP
                Re(1,var+1)              = Ri(1,var)-INS_grades;   % Raggio esterno grade (i+1)-esimo del WP     
                WP_w0(1,var)             = Cond_w(1,var)*n_turns(1,var); % massimo ingombro toroidale WP
                check_w = 2*Ri(1,var)*tan(theta_TF/2);
                while (check_w-(WP_w0(1,var)+GoundIns*2))/2 <= toroidal_gap
                        n_turns(1,var) = n_turns(1,var)-2;     
                        % Definisco dimensioni interne del cavo CICC per l'i-esimo grade
                        r_SC(1,var)              = (0.005)*Increm; %JT(1,var);   % Impogno raggio curvatura corner cavo = a spessore Jakcet
                        tins(1,var)               = (0.001)*Increm; % Isolante di spira
                        Cond_w(1,var)            = WP_w0(1,1)/n_turns(1,1);  
                    if cell2sym(type_cable(1,var)) == 'HTS'  
                        E_cbl = E_cbl_HTS;
                    else
                        E_cbl = E_cbl_LTS;
                    end	
                    if shape_cable == 200
                        SC_w(1,var)              = 2*sqrt(S_Cable(1,var)/pi);
                        SC_h(1,var)              = SC_w(1,var);   % Ottengo altezza cavo SC
                        R_J(1,var)               = r_SC(1,var);   % Raggio di curvatura corner Jacket
                        Cond_h(1,var)            = Cond_w(1,var); % Ottengo l'altezza del cavo RIS_ 
                        JT(1,var)                = (Cond_w(1,var)-2*tins(1,var)-SC_w(1,var))/2;
                        Ke_cavo_rad(1,var) = 2*E_jckt*JT(1,var)/Cond_h(1,var)+2*tins(1,var)*E_ins/Cond_h(1,var)+...
                        +(1/(E_cbl*SC_w(1,var)/SC_h(1,var))+2/(E_jckt*Cond_w(1,var)/JT(1,var))+2/(E_ins*Cond_w(1,var)/tins(1,var)))^-1;
                        Ke_cavo_tor(1,var) = 2*E_jckt*JT(1,var)/Cond_w(1,var)+2*tins(1,var)*E_ins/Cond_w(1,var)+...
                        +(1/(E_cbl*SC_h(1,var)/SC_w(1,var))+2/(E_jckt*JT(1,var)/Cond_w(1,var))+2/(E_ins*tins(1,var)/Cond_w(1,var)))^-1;
                    elseif shape_cable == 201
                        JT(1,var) =  min_JT;
                        S_rm_grades(1,var) = 1e30;
                    while S_rm_grades(1,var) > S_amm_JT/1.3
                        JT(1,var) = JT(1,var) + 1e-4;
                        SC_w(1,var)              = Cond_w(1,var)-JT(1,var)*2-tins(1,var)*2; % Ottengo larghezza cavo SC 
                        SC_h(1,var)              = (S_Cable(1,var)+(4-pi)*r_SC(1,var)^2)/SC_w(1,var); % Ottengo altezza cavo SC
                        R_J(1,var)               = r_SC(1,var) + JT(1,var); % Raggio di curvatura corner Jacket
                        Cond_h(1,var)            = SC_h(1,var)+JT(1,var)*2+tins(1,var)*2; % Ottengo l'altezza del cavo Rect                        
                        Ke_cavo_rad(1,var) = 2*E_jckt*JT(1,var)/Cond_h(1,var)+2*tins(1,var)*E_ins/Cond_h(1,var)+...
                        +(1/(E_cbl*SC_w(1,var)/SC_h(1,var))+2/(E_jckt*Cond_w(1,var)/JT(1,var))+2/(E_ins*Cond_w(1,var)/tins(1,var)))^-1;
                        Ke_cavo_tor(1,var) = 2*E_jckt*JT(1,var)/Cond_w(1,var)+2*tins(1,var)*E_ins/Cond_w(1,var)+...
                        +(1/(E_cbl*SC_h(1,var)/SC_w(1,var))+2/(E_jckt*JT(1,var)/Cond_w(1,var))+2/(E_ins*tins(1,var)/Cond_w(1,var)))^-1;
                        K_jckt = 2*JT(1,var)/Cond_h(1,var)*E_jckt;
                        dcr_jckt = K_jckt/Ke_cavo_rad(1,var);
                        r_steel = (Cond_w(1,var)-2*tins(1,var))/(2*JT(1,var));
                        S_rm_grades(1,var) = p_rs*r_steel*dcr_jckt+S_z_JT;      % Radial memebrane stress Jacket innermost layer 
                    end  
                    end
                    S_CICC(1,var)            = ((Cond_w(1,var)-2*tins(1,var))*(Cond_h(1,var)-2*tins(1,var)))-((4-pi)*R_J(1,var)^2); % Sezione cavo non isolato
                    S_JT(1,var)              = S_CICC(1,var)-S_Cable(1,var); % Sezione di acciaio
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
                    S_Cable(1,var)           = S_Cable(1,var-1);
                    type_cable(1,var)        = type_cable(1,var-1);
                    % Definisco dimensioni interne del cavo CICC per l'i-esimo grade
                    r_SC(1,var)              = (0.005)*Increm; %JT(1,var);   % Impogno raggio curvatura corner cavo = a spessore Jakcet
                    tins(1,var)               = (0.001)*Increm; % Isolante di spira
                    Cond_w(1,var)            = WP_w0(1,1)/n_turns(1,1);  
                    if cell2sym(type_cable(1,var)) == 'HTS'  
						E_cbl = E_cbl_HTS;
					else
						E_cbl = E_cbl_LTS;
                    end	
                    if shape_cable == 200
                        SC_w(1,var)              = 2*sqrt(S_Cable(1,var)/pi);
                        SC_h(1,var)              = SC_w(1,var);   % Ottengo altezza cavo SC
                        R_J(1,var)               = r_SC(1,var);   % Raggio di curvatura corner Jacket
                        Cond_h(1,var)            = Cond_w(1,var); % Ottengo l'altezza del cavo RIS_ 
                        JT(1,var)                = (Cond_w(1,var)-2*tins(1,var)-SC_w(1,var))/2;
                        Ke_cavo_rad(1,var) = 2*E_jckt*JT(1,var)/Cond_h(1,var)+2*tins(1,var)*E_ins/Cond_h(1,var)+...
                                     +(1/(E_cbl*SC_w(1,var)/SC_h(1,var))+2/(E_jckt*Cond_w(1,var)/JT(1,var))+2/(E_ins*Cond_w(1,var)/tins(1,var)))^-1;
                        Ke_cavo_tor(1,var) = 2*E_jckt*JT(1,var)/Cond_w(1,var)+2*tins(1,var)*E_ins/Cond_w(1,var)+...
                                     +(1/(E_cbl*SC_h(1,var)/SC_w(1,var))+2/(E_jckt*JT(1,var)/Cond_w(1,var))+2/(E_ins*tins(1,var)/Cond_w(1,var)))^-1;
                    elseif shape_cable == 201 
                        JT(1,var) = min_JT;
                        S_rm_grades(1,var) = 1e30;
                        while S_rm_grades(1,var) > S_amm_JT/1.3
                            JT(1,var) = JT(1,var) + 1e-4;
                            SC_w(1,var)              = Cond_w(1,var)-JT(1,var)*2-tins(1,var)*2; % Ottengo larghezza cavo SC 
                            SC_h(1,var)              = (S_Cable(1,var)+(4-pi)*r_SC(1,var)^2)/SC_w(1,var); % Ottengo altezza cavo SC
                            R_J(1,var)               = r_SC(1,var) + JT(1,var); % Raggio di curvatura corner Jacket
                            Cond_h(1,var)            = SC_h(1,var)+JT(1,var)*2+tins(1,var)*2; % Ottengo l'altezza del cavo Rect                        
                            Ke_cavo_rad(1,var) = 2*E_jckt*JT(1,var)/Cond_h(1,var)+2*tins(1,var)*E_ins/Cond_h(1,var)+...
                                     +(1/(E_cbl*SC_w(1,var)/SC_h(1,var))+2/(E_jckt*Cond_w(1,var)/JT(1,var))+2/(E_ins*Cond_w(1,var)/tins(1,var)))^-1;
                            Ke_cavo_tor(1,var) = 2*E_jckt*JT(1,var)/Cond_w(1,var)+2*tins(1,var)*E_ins/Cond_w(1,var)+...
                                     +(1/(E_cbl*SC_h(1,var)/SC_w(1,var))+2/(E_jckt*JT(1,var)/Cond_w(1,var))+2/(E_ins*tins(1,var)/Cond_w(1,var)))^-1;
                            K_jckt = 2*JT(1,var)/Cond_h(1,var)*E_jckt;
                            dcr_jckt = K_jckt/Ke_cavo_rad(1,var);
                            r_steel = (Cond_w(1,var)-2*tins(1,var))/(2*JT(1,var));
                            S_rm_grades(1,var) = p_rs*r_steel*dcr_jckt+S_z_JT;      % Radial memebrane stress Jacket innermost layer 
                        end  
                    end
                    S_CICC(1,var)            = ((Cond_w(1,var)-2*tins(1,var))*(Cond_h(1,var)-2*tins(1,var)))-((4-pi)*R_J(1,var)^2); % Sezione cavo non isolato
                    S_JT(1,var)              = S_CICC(1,var)-S_Cable(1,var); % Sezione di acciaio                                
                    Ri(1,var)                = Re(1,var)-Cond_h(1,var); % Raggio interno grade i-esimo del WP
                    Re(1,var+1)              = Ri(1,var)-INS_grades;    % Raggio esterno grade (i+1)-esimo del WP 
                    WP_w0(1,var)             = Cond_w(1,var)*n_turns(1,var); % massimo ingombro toroidale WP
                end  
            end  
            if n_layers > maxdim || n_layers < 0
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
            WP_w = WP_w0(1,1);
            WP_w_tot = WP_w+2*GoundIns;
            WP_h = sum(Cond_h(1:n_layers));
            WP_h_tot = WP_h+2*GoundIns;
            A_WP = sum(Cond_h(1:n_layers).*Cond_w(1:n_layers).*n_turns(1:n_layers));
            A_CICC_tot = sum(S_CICC(1:n_layers).*n_turns(1:n_layers));  % Sezione totale cavo non isolato
            A_SC_tot = sum(S_Cable(1:n_layers).*n_turns(1:n_layers));   % Sezione totale S/C+stabilizer 
            A_JT_tot = sum(S_JT(1:n_layers).*n_turns(1:n_layers));      % Sezione totale acciaio nei jacket  
            Ri_ = R_TF_Innerleg;                            % Quota a raggio esterno Case
            Rj_ = Ri_-WP_h-dr_plasma_side-GoundIns*2;       % Quota a raggio interno WP con Ground     

            check_w(1,1:n_layers) = 2*Ri(1,1:n_layers).*tan(theta_TF/2);          % massimo ingombro toroidale Case
            if min((check_w-(WP_w0(1,1:n_layers)+GoundIns*2))/2) < toroidal_gap
                continue
            end

%% Cechck geometrico sulle dimensioni ottenute nel cavo
            r_cable(1:n_layers) = Cond_w(1:n_layers)./Cond_h(1:n_layers);             % Aspect ratio cavi di ogni grades
            if min(SC_w) <= 0.005 || min(r_cable(1:n_layers))< 0.99 || min(JT(1:n_layers)) < min_JT
                continue
            end

%% Primary radial stress (Pm+Pb)
            p_rs = B_TF^2/(2*Mu_0);                                 % Magnetic pressure, thin WP         
            param = Cond_w(1,n_layers)/SC_w(1,n_layers);
            if param > 1 && param < 2.8
                xxx = [1.087,1.136,1.190,1.250,1.316,1.389,1.471,1.563,1.667,1.786,1.923,2.083,2.273,2.500,2.778]; 
                if cell2sym(type_cable(1,var)) == 'LTS'
                    E_cbl = E_cbl_LTS;
                    yyy = [1.01,1.03,1.06,1.10,1.16,1.22,1.29,1.36,1.43,1.50,1.57,1.64,1.71,1.78,1.85];  
                elseif cell2sym(type_cable(1,var)) == 'HTS'
                    E_cbl = E_cbl_HTS;
                    yyy = [1.01,1.02,1.02,1.04,1.06,1.07,1.11,1.13,1.16,1.18,1.22,1.27,1.29,1.29,1.31]; 
                end
                p = polyfit(xxx,yyy,5);
                scf = polyval(p,param);        
            else
                scf = 1.5;
            end
            var = n_layers;
            Ke_cavo_rad(1,var) = 2*E_jckt*JT(1,var)/Cond_h(1,var)+2*tins(1,var)*E_ins/Cond_h(1,var)+...
                                 +(1/(E_cbl*SC_w(1,var)/SC_h(1,var))+2/(E_jckt*Cond_w(1,var)/JT(1,var))+2/(E_ins*Cond_w(1,var)/tins(1,var)))^-1;
            Ke_cavo_tor(1,var) = 2*E_jckt*JT(1,var)/Cond_w(1,var)+2*tins(1,var)*E_ins/Cond_w(1,var)+...
                                 +(1/(E_cbl*SC_h(1,var)/SC_w(1,var))+2/(E_jckt*JT(1,var)/Cond_w(1,var))+2/(E_ins*tins(1,var)/Cond_w(1,var)))^-1;
            K_jckt = 2*JT(1,var)/Cond_h(1,var)*E_jckt;
            dcr_jckt = K_jckt/Ke_cavo_rad(1,var);
            r_steel = (Cond_w(n_layers)-2*tins(n_layers))/(2*JT(n_layers));
            S_rm = p_rs*r_steel*scf*dcr_jckt; % Radial memebrane stress Jacket innermost layer

%% 
            clear Rk_ DTF S_T_JT S_T_VT Ke_WP_rad Ke_WP_tor K_ps_rad K_ps_tor dcr_vault_tor dcr_WP_rad K_vault_rad K_vault_tor K_lat_rad K_lat_tor
            DTF = 0.05;    % Vaul width 
            S_T_JT = 1e30;
            S_T_VT = 1e30; % Test value for vault Tresca stress
            while S_T_VT > S_amm_VT || S_c_VT > S_amm_VT/1.3 || S_T_JT > S_amm_JT || S_rm_JT > S_amm_JT/1.3
                DTF = DTF+0.001; % Se supero Tresca incremento spessore naso TF
                Rk_ = Rj_-DTF; % Innermost Case radius 
                CASE_w_l = 2*Rk_*tan(theta_TF/2); % Case low part width 
                A_tot = (CASE_w+CASE_w_l)*(Ri_-Rk_)/2;
                A_CASE = A_tot-A_WP;
                A_VT = (2*Rj_*tan(theta_TF/2)+CASE_w_l)*(Rj_-Rk_)/2; % Vault section trp
                A_VT_circ = pi*(Rj_^2-Rk_^2)*1/n_TF; % Vault section circ
                % Rigidezze radiali 
                Ke_WP_rad = sum(1./(Ke_cavo_rad(1:n_layers).*n_turns(1:n_layers)))^-1;        
                K_ps_rad = E_case/dr_plasma_side*(2*Ri_*tan(theta_TF/2));
                K_lat_rad = E_case*(lateral_w/2)/WP_h;
                K_vault_rad = (E_case*CASE_w_l/DTF);   
                Ke_case_rad = (1/(Ke_WP_rad+2*K_lat_rad)+1/K_vault_rad+1/K_ps_rad)^-1;
                dcr_WP_rad = Ke_WP_rad/Ke_case_rad; 
                S_rm_JT = S_rm*dcr_WP_rad; % Radial memebrane stress Jacket innermost layer correction          
                
                h_unit = 1;
                k_steel_tor = E_case*(h_unit*(Ri_-Rk_)*2*pi*(Ri_+Rk_)/2); % full steel casing
                k_SC_tor = E_cbl*(h_unit*(Ri_-Rj_)*2*pi*(Ri_+Rj_)/2)*(1-(A_WP-A_SC_tot)/A_WP);  
                k_JT_tor = E_jckt*(h_unit*(Ri_-Rj_)*2*pi*(Ri_+Rj_)/2)*(1-(A_WP-A_JT_tot)/A_WP); 
                k_vault_tor = E_case*(h_unit*(Rj_-Rk_)*2*pi*(Rj_+Rk_)/2);
            
                dcr_vault_tor = k_steel_tor/(k_SC_tor+k_JT_tor+k_vault_tor);
            
                beta = Rk_/Ri_;
                S_c_VT = 2/(1-beta^2)*p_rs*dcr_vault_tor; % @ Rk_ correggere con 
            
                k_bf = 0.5*log(RTFo(dp)/RTFi(dp)); % k bending free
                T_bf = 0.5*(k_bf*n_TF*(n_spire_(1,1)*Iop)^2*Mu_0/(2*pi)); % Hoop tension along TF longitudinal axis
            
                S_z =  T_bf/(A_JT_tot+A_CASE);  
            
                S_T_VT = (S_z+S_c_VT); % Vault Tresca stress 
                S_T_JT = (S_z+S_rm_JT);  % Jacket Tresca stress      
            end   
            V = L*Iop/Tau_discharge*1e-3; % Tensione max singolo TF
%%
            if S_T_VT < S_amm_JT && S_T_JT < S_amm_JT && S_T_JT>0 && S_T_VT>0 
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
                VV_BB = R_0-Ri-R_0/A(dp);
                DATA(counter,:) = table(dp,S_T_VT,S_T_JT,R_0,B_PHI_0,B_TF,Iop,JENG,L,E,Ri_,Rj_,Rk_,radial_build,Nose,WP_h,WP_w,...
                    lateral_w,n_cond,n_layers,n_turns,type_cable,Cond_w,Cond_h,JT,r_cable,N_Sc,N_Cu,S_Cable,S_REBCO,S_Cu_HTS,THS,Tau_discharge);
                plot_data(counter,:) = [S_T_VT,S_T_JT,Rk_,Iop,B_PHI_0,B_TF];                                                                                                    
            end             
        end        
    end  
end

%% Scrivo matrice soluzione prodotte
if size(DATA,1) ~=0     
   writetable(DATA,TitleName) 
end

%%
figure('units','normalized','outerposition',[0 0 1 1])
Ii = table2array((DATA(1:end,7)));
Bb = table2array((DATA(1:end,6)));
Rkk = table2array((DATA(1:end,13)));
Sint = table2array((DATA(1:end,2)));
Je = table2array((DATA(1:end,8)));
scatter(Ii.*1e-3,Rkk,20,Sint,'filled') 
colormap(jet);
set(gca,'FontSize',20);
colorbar
c = colorbar;
c.Label.FontSize = 20;
c.Label.Color = 'k';
c.Label.Rotation = 90;
c.Label.String = 'SINT JT [MPa]';
c.Label.Interpreter = 'latex';
c.TickLabelInterpreter = 'latex';
mod_max = round(max(Sint));
mod_min = round(min(Sint));
clim([mod_min mod_max]);
MAXMIN=get(c,'Limits');
T = (linspace(MAXMIN(1),MAXMIN(2),9));
set(c,'Ticks',T)
ylabel('Rk [m]','fontsize',20,'Interpreter','latex')
xlabel('Iop [kA]','fontsize',20,'Interpreter','latex')
TitleName = sprintf('Iop vs Rk - A%G',A(dp)*10);
title(TitleName,'fontsize',20,'Interpreter','latex')
grid on
print('-dpng','-r300',TitleName);
%%
figure('units','normalized','outerposition',[0 0 1 1])
scatter(Je,Rkk,20,Sint,'filled') 
[x1,S] = polyfit(Je,Rkk,2);
[y1,delta] =  polyval(x1,min(Je):1e-3:max(Je),S);
hold on
plot(min(Je):1e-3:max(Je),y1,'k','linewidth',1)

colormap(jet);
set(gca,'FontSize',20);
colorbar
c = colorbar;
c.Label.FontSize = 20;
c.Label.Color = 'k';
c.Label.Rotation = 90;
c.Label.String = 'SINT JT [MPa]';
c.Label.Interpreter = 'latex';
c.TickLabelInterpreter = 'latex';
mod_max = round(max(Sint));
mod_min = round(min(Sint));
clim([mod_min mod_max]);
MAXMIN=get(c,'Limits');
T = (linspace(MAXMIN(1),MAXMIN(2),9));
set(c,'Ticks',T)
ylabel('Rk [m]','fontsize',20,'Interpreter','latex')
xlabel('Jeng [$A/mm^2$]','fontsize',20,'Interpreter','latex')
TitleName = sprintf('Rk vs Jeng - A%G',A(dp)*10);
title(TitleName,'fontsize',20,'Interpreter','latex')
grid on
print('-dpng','-r300',TitleName);
%% 
% figure('units','normalized','outerposition',[0.25 0.1 0.5 0.9]) Ii = table2array((DATA(1:end,6))); 
% U = table2array((DATA(1:end,31))); Bb = table2array((DATA(1:end,5))); Rkk = 
% table2array((DATA(1:end,9))); Sint = table2array((DATA(1:end,3))); NIi = table2array((DATA(1:end,6))).*table2array((DATA(1:end,15)))*1e-6;
% 
% scatter(U,Rkk,20,Sint,'filled') colormap(jet); set(gca,'FontSize',20); colorbar 
% c = colorbar; c.Label.FontSize = 20; c.Label.Color = 'k'; c.Label.Rotation = 
% 90; c.Label.String = 'SINT JT [MPa]'; c.Label.Interpreter = 'latex'; c.TickLabelInterpreter 
% = 'latex'; mod_max = round(max(Sint)); mod_min = round(min(Sint)); caxis([mod_min 
% mod_max]); MAXMIN=get(c,'Limits'); T = (linspace(MAXMIN(1),MAXMIN(2),8)); set(c,'Ticks',T) 
% ylabel('Rk [m]','fontsize',20,'Interpreter','latex') xlabel('U [$MJ$]','fontsize',20,'Interpreter','latex') 
% TitleName = sprintf('Rk vs U - A%G WP %s',A(dp)*10,WP_SC_type); title(TitleName,'fontsize',20,'Interpreter','latex') 
% grid on print('-dpng','-r300',TitleName);

% % % load handel
% % % sound(y,Fs)