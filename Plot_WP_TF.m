clearvars; clear memory;  clc; close all;

%%
printcond = 0;
%%
Increm = 1.000; % Per passare a dimensioni a Tamb (1.003)
counter = 0;
iter = 0;
%% Design point input
n_layers = 9;
maxdim = n_layers;
n_turns = ones(1,n_layers)*[12];
% n_turns = [10 10 10 10 10 10 10 10 10 10 ];
n_grades = 1;
grading_at = [1 (n_layers+1)]; %% il grade viene applicato a partire dal layer indicato (incluso)
%%
lateral_w = 0.07; 
min_JT = 0.002; 
V_MAX = 5e3/2; % [V]
%%
dp = 1;
n_TF = 12; 
%% input da esplorazione 
B0(dp) = 5.4;  	     
R0(dp) = 2.53;  
A(dp) = 4.6;
wb = 0.80;
RTFi(dp) = R0(dp)-R0(dp)/A(dp)-wb;    % [m]
ripple = 0.006;
RTFo(dp) = (R0(dp)+R0(dp)/A(dp))*(1/ripple)^(1/n_TF);    % [m]
R_VV = RTFi(dp)*1.05;  
%%
corr_B_WP   = 1.08;
WP_SC_type  = 'full_HTS'; % full_LTS - full_HTS - Hybrid__
shape_cable = 'Rect'; % RIS_ - Rect
%% input da design points
% R0   = [7.5 6.5 7.4 9.0 3.2 8.9380 2.55 2.5 7.8 7.2 7.4 8.3];           % [m]
% B0   = [4.0 6.6 11.4 5.0 7.4 4.89 6.0 7.0 3.9 5.9 8.1 13.25];               % [T]
% A    = [2.6 3.3 4.5 3.1 6.0 3.1 5.3 5.5 2.6 3.1 4.0 4.6];
% RTFi = [3.02 2.9 4.178 4.3 1.77 4.6248 1.22 1.295 3.4 3.57 4.05 5.29];       % [m]
% RTFo = [13.837 11.9 12.310 12.0 5.26  16.2756 4.456 4.39 13.25 11.87 12.20 12.20];  % [m]
% H_TF = [19.9 16.4 15 0 0];  
%% Ammissibili acciaio
S_amm_JT                   = (1000*1e6);
S_amm_VT                   = (1000*1e6);
S_amm_Cu                   = (100*1e6);
%% Parametri operativi
Mu_0 = 4e-7*pi;
theta_TF = 2*pi/n_TF;
amp_corr_R0 = 1/(R0(dp)/RTFi(dp))*(1+1/((R0(dp)/RTFi(dp))^n_TF-1)+1/((RTFo(dp)/R0(dp))^n_TF-1));
NI = (2*pi*R0(dp)*B0(dp)/Mu_0)/n_TF*1e-6;        % Total TF current [MA]
B_PHI_TF = Mu_0*n_TF*NI*1e6/(2*pi*RTFi(dp));     % Max field on TF
B_PHI_0 = B_PHI_TF*amp_corr_R0;
B_PHI_TF = B_PHI_TF*corr_B_WP;   % Max field on TF correction
NI = NI*1e6;                              % Total TF current [A]
R_TF_Outerleg = RTFo(dp);                 % Outer-leg inner radius
R_TF_Innerleg = RTFi(dp);                 % Inner-leg outer radius
S_VV = 90*1e6;                           % VV Yeld limit
%%
Tau_discharge1 = B0(dp)*NI*n_TF*(R0(dp)/A(dp))^2/(R_VV*S_VV);
%% materials
E_jckt = 205;
E_cbl_HTS = 120;
E_ins = 12;
E_cbl_LTS = 0;
E_case = 205;
%% WP dimensioning
CASE_w          = 2*R_TF_Innerleg*tan(theta_TF/2); % Case width
dr_plasma_side  = 0.02; % [m] Spessore case fronte plasma
GoundIns        = 0.005; % Ground insulation WP
WP_w            = 2*(R_TF_Innerleg-dr_plasma_side-GoundIns)*tan(theta_TF/2)-2*lateral_w; % WP width
INS_grades      = 0.0005; % Spessore isolante tra i Grades
%% Dati cavo superconduttore Nb3Sn
d_fili = 0.82; % wire diameter [mm] 
cos_theta = 0.97; % Cos(theta): considers that wires are twisted, then the real effective section is grater 
VF = 0.8 ;       % Void Fraction CICC 27%
r_cc = 0.005/2;     % Cooling channel radius
N_fili = [1500 1440 1350 1296 1200 1152 1080 972 960 900 864 810 768 720 675 648 540 486 420 360 324 300 216 180 162 144]; % Possible triplets combinations
%% Dati configurazioni da iterare 
iter = iter+1;                                       
%% Turns in each grade
% n_layers_ = ones(1,n_grades);
n_spire_ = zeros(1,n_layers);
n_spire_(1,1) = sum(n_turns);
for var = 2:n_layers 
    n_spire_(1,var) = n_spire_(1,var-1)-n_turns(1,var-1);
end   
%% Preallocations
Cond_h = zeros(1,maxdim); 
Cond_w = zeros(1,maxdim); 
clear EQV_Cable_Area;
% EQV_Cable_Area = zeros(1,maxdim); 
JT = zeros(1,maxdim); 
r_cable = zeros(1,maxdim); 
for k = 1:maxdim
    type_cable(1,k) = {'---'};
end 
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
%% Definition of the magnetic field peak in each grade (linera beahvior from B(Re)= Bmax to B(Ri) = 0                 
Re(1,1) = R_TF_Innerleg-dr_plasma_side-GoundIns; % WP innerl-leg outer radius (non-insulated)
B_TF = B_PHI_TF; % (n_TF*n_spire_(1,1)*Iop*Mu_0)/(2*pi*Re(1,1)); 
B_layers = B_TF.*(n_spire_./n_spire_(1,1)); % Ottengo B considerando un andamento lineare nel WP
% xB = [10 12 14 16.6];
% yI = [1420 800 580 200];
% p = polyfit(xB,yI,2);
clear Ic_sc_HJc Ic_sc
for var = 1:n_layers                                                                                
    B = B_layers(1,var);
    Ic_sc(1,var) = Ic_Nb3Sn_WST(B,d_fili); % Critical current Nb3Sn
    % Ic_sc(1,var) = Ic_NbTi_TF(B,d_fili); % Critical current Nb3Sn  
    % Ic_sc_HJc(1,var) = polyval(p,B);
end     
N_sc0 = Iop./Ic_sc; % Needed wires in LTS cable for each grade              
%% solo HTS
if WP_SC_type == 'full_HTS'
   LTS_Bmin = 0;
elseif WP_SC_type == 'full_LTS' 
   LTS_Bmin = 15.5;
elseif WP_SC_type == 'Hybrid__' 
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
% if WP_SC_type == 'full_LTS'
%     if  min(check_type_HTS(check_type_HTS>0)) >0
%         continue
%     end
% end
%% Induttanza TF modello shell
clear L Tau_discharge
% L = Mu_0*R0(dp)*(n_TF*n_spire_(1,1))^2*(1-sqrt(1-(R0(dp)-RTFi(dp))/R0(dp)))/n_TF*1.1;
% 
k = 0.5*log(RTFo(dp)/RTFi(dp));
r0 = sqrt(RTFo(dp)*RTFi(dp));
L = (1/2)*(Mu_0*r0*(n_TF*n_spire_(1,1))^2*(k^2)*(besseli(0,k)+2*besseli(1,k)+besseli(2,k)))/n_TF;
%
Tau_discharge2 = (L*Iop/V_MAX); % [s] - scarico in gruppi di tre bobine
Tau_discharge = max([Tau_discharge1 Tau_discharge2 4]);
E = 1/2*L*n_TF*Iop^2*1e-6;
%% 
for n = 1:n_grades
    var = grading_at(n);
%     var_LTS = min(check_type_LTS(check_type_LTS>0)):max(check_type_LTS); 
%     var = min(check_type_LTS(check_type_LTS>0));
%     if min(check_type_LTS(check_type_LTS>0)) >0 
    if check_type_LTS(var) >0 
        % LTS
        E_cbl = E_cbl_LTS;                 
        a_ = N_fili > N_sc0(1,var);
        z_ = N_fili(a_);
        N_sc(1,var) = min(z_);
        % N_sc(1,var) = N_sc0(1,var);
        S_Cu(1,var) = N_sc(1,var)*pi*d_fili^2/8; % Sezione di rame non segregato 
        Cu_cross_sect_in_SC(1,var) = S_Cu(1,var)*1e-6; % Sezione di rame non segregato in [m2]
        Nb3Sn_cross_section(1,var) = S_Cu(1,var)*1e-6;  % controllare rapporto Cu non Cu; qui si considera 1
        J_copper = 180; % Densità di corrente critica stabilizzatore [A/mm2]
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
        THS(1,var) = Ths_TF_Nb3Sn(Iop,Bp,Cus,Cuin,NB3SN,Tau_discharge);    
        while THS(1,var)>250 && N_Cu(1,var)>0 
            N_Cu(1,var) = N_Cu(1,var)+1;
            Bp = B_layers(1,var);                           
            Cus = pi*N_Cu(1,var)*d_fili^2/4*1e-6;  
            Cuin = Cu_cross_sect_in_SC(1,var);
            NB3SN = Nb3Sn_cross_section(1,var);
            THS(1,var) = Ths_TF_Nb3Sn(Iop,Bp,Cus,Cuin,NB3SN,Tau_discharge);                     
        end
        while THS(1,var)<250 && N_Cu(1,var)>0 && THS(1,var) > 0
            N_Cu(1,var) = N_Cu(1,var)-1;
            Bp = B_layers(1,var);                           
            Cus = pi*N_Cu(1,var)*d_fili^2/4*1e-6;  
            Cuin = Cu_cross_sect_in_SC(1,var);
            NB3SN = Nb3Sn_cross_section(1,var);
            THS(1,var) = Ths_TF_Nb3Sn(Iop,Bp,Cus,Cuin,NB3SN,Tau_discharge);                      
        end
        N_tot(1,var) = N_sc(1,var)+N_Cu(1,var);
        EQV_Cable_Area_0       = (N_tot(1,var)*pi*d_fili^2/4/cos_theta/VF+(pi*r_cc^2))*1e-6; % Area eqv cavo LTS
        D_eqv                  = (EQV_Cable_Area_0*4/pi)^0.5; % Diametro cavo equivalente da area eqv
        A_w                    = ((D_eqv+0.0004)^2-D_eqv^2)*pi/4; % Considereo l'area del wrapping in acciaio di spessore 0.4 mm
        EQV_Cable_Area(1,var)  =  pi/4*D_eqv^2+A_w; % Correggo area equivalente del cavo LTS
        %
    %     THS(1,var_LTS)    = THS(1,var);
    %     N_Cu(1,var_LTS)  = N_Cu(1,var);
    %     N_sc(1,var_LTS)    = N_sc(1,var);
    %     EQV_Cable_Area(1,var_LTS)    = EQV_Cable_Area(1,var);
        step_ = grading_at(n):grading_at(n+1)-1;
        THS(1,step_)    = THS(1,var);
        N_Cu(1,step_)  = N_Cu(1,var);
        N_sc(1,step_)    = N_sc(1,var);
        EQV_Cable_Area(1,step_)    = EQV_Cable_Area(1,var);
    end        
% var_HTS = min(check_type_HTS(check_type_HTS>0)):max(check_type_HTS); 
% var = min(check_type_HTS(check_type_HTS>0));
% if  min(check_type_HTS(check_type_HTS>0)) >0  
    if  check_type_HTS(var) >0  
        %% HTS                     
        E_cbl = E_cbl_HTS;
        type_cable(1,var)        = {'HTS'};
        % % % JB_REBCO                = (0.9000*(B_layers(1,var))^2-54.100*(B_layers(1,var))+1.1455e3)*1e6; % Densità di corrente critica REBCO
        JB_REBCO                 = Jc_REBCO(B_layers(1,var),S_tapes)*1e6;
        S_REBCO(1,var)           = Iop/JB_REBCO;
        J_copper                 = 120e6; % Densità di corrente critica stabilizzatore [A/m2]
        S_Cu_HTS(1,var)          = Iop/J_copper;        
        Bp = B_layers(1,var);                           
        Cus = S_Cu_HTS(1,var); % Sezione di rame segregato 
        Cuin = 0; % Sezione di rame non segregato
        NB3SN = 0;%S_REBCO(1,var);
        THS(1,var) = Ths_TF_Nb3Sn(Iop,Bp,Cus,Cuin,NB3SN,Tau_discharge);   
        while THS(1,var)>250
            S_Cu_HTS(1,var) = S_Cu_HTS(1,var)*1.01;
            Bp = B_layers(1,var);                           
            Cus = S_Cu_HTS(1,var);  
            Cuin = 0;
            NB3SN = 0;%S_REBCO(1,var);
            THS(1,var) = Ths_TF_Nb3Sn(Iop,Bp,Cus,Cuin,NB3SN,Tau_discharge);                      
        end                  
        while THS(1,var)<250 && S_Cu_HTS(1,var)>0
            S_Cu_HTS(1,var) = S_Cu_HTS(1,var)*0.999;
            Bp = B_layers(1,var);                           
            Cus = S_Cu_HTS(1,var);  
            Cuin = 0;
            NB3SN = 0;%S_REBCO(1,var);
            THS(1,var) = Ths_TF_Nb3Sn(Iop,Bp,Cus,Cuin,NB3SN,Tau_discharge);                       
        end                                         
        EQV_Cable_HTS_0          = S_REBCO(1,var)+S_Cu_HTS(1,var)+(pi*r_cc^2); % Area eqv cavo HTS 
        D_eqv(1,var)             = (EQV_Cable_HTS_0*4/pi)^0.5; % Diametro cavo equivalente da area eqv
        A_wr(1,var)              = ((D_eqv(1,var)+0.0005)^2-D_eqv(1,var)^2)*pi/4;% Jacket tondo interno (design cavo HTS ENEA)
        EQV_Cable_Area(1,var)    = EQV_Cable_HTS_0+A_wr(1,var); % Correggo area equivalente del cavo HTS
        %
        % % % EQV_Cable_Area(1,var)*1e6
        % % % S_REBCO(1,var)*1e6
        % % % S_Cu_HTS(1,var)*1e6
%         THS(1,var_HTS)    = THS(1,var);
%         S_Cu_HTS(1,var_HTS)    = S_Cu_HTS(1,var);
%         S_REBCO(1,var_HTS)    = S_REBCO(1,var);
%         EQV_Cable_Area(1,var_HTS)    = EQV_Cable_Area(1,var); 
        step_ = grading_at(n):grading_at(n+1)-1;
        THS(1,step_)    = THS(1,var);
        N_Cu(1,step_)  = N_Cu(1,var);
        N_sc(1,step_)    = N_sc(1,var);
        EQV_Cable_Area(1,step_)    = EQV_Cable_Area(1,var);
    end  
end
%%  Definisco sezioni cavi in ogni layers 
R_WP_IL = R_TF_Innerleg;
R_WP_OL = R_TF_Outerleg;
k_bf = 0.5*log(R_WP_OL/R_WP_IL); % k bending free
T_bf = 0.5*(k_bf*n_TF*(n_spire_(1,1)*Iop)^2*Mu_0/(2*pi)); % Hoop tension along TF longitudinal axis
T_bf_  = T_bf*1e-6;
Fr = (n_TF*(n_spire_(1,1)*Iop)^2*Mu_0/2)*(1-(1/sqrt(1-(RTFi(dp)/R0(dp))^2)));
F0 = n_TF*(n_spire_(1,1)*Iop)^2*Mu_0/4/pi;
H = 1;
Fr = -F0*H/RTFi(dp)*1e-6;

WP_w0(1,1) = 2*Re(1,1)*tan(theta_TF/2)-lateral_w*2-GoundIns*2; % massimo ingombro toroidale WP
S_z_JT =  T_bf/(WP_w0(1,1)^2)/2;
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
        JT(1,var) = min_JT;
        S_rm_grades(1,var) = 1e30;
        while S_rm_grades(1,var) > S_amm_JT/1.3
            JT(1,var) = JT(1,var)+5e-4;
            SC_w(1,var)              = Cond_w(1,var)-JT(1,var)*2-INS(1,var)*2; % Ottengo larghezza cavo SC 
            SC_h(1,var)              = (EQV_Cable_Area(1,var)+(4-pi)*R_SC(1,var)^2)/SC_w(1,var); % Ottengo altezza cavo SC
            R_J(1,var)               = R_SC(1,var) + JT(1,var); % Raggio di curvatura corner Jacket
            Cond_h(1,var)            = SC_h(1,var) + JT(1,var)*2+INS(1,var)*2; % Ottengo l'altezza del cavo Rect                        
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
end
%% Cechck geometrico sulle dimensioni ottenute nel cavo
r_cable(1:n_layers) = Cond_w(1:n_layers)./Cond_h(1:n_layers);             % Aspect ratio cavi di ogni grades
%% Dati WP    
WP_w = WP_w0(1,1);
WP_h = Re(1,1)-Ri(1,n_layers);                % Altezza del WP
A_WP = sum(Cond_h(1:n_layers).*WP_w0(1:n_layers));                     % Area WP considerando un rettangolo senza corner 
A_SC_tot = A_SC(1)*n_spire_(1);               % Sezione totale superconduttore 
A_J_tot = A_J(1)*n_spire_(1);                 % Sezione totale acciaio nei jacket  
Ri_ = R_TF_Innerleg;                          % Quota a raggio esterno Case
Rj_ = Ri_-WP_h-dr_plasma_side-GoundIns*2;     % Quota a raggio interno WP con Ground                              
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
Ke_cavo_rad(1,var) = 2*E_jckt*JT(1,var)/Cond_h(1,var)+2*INS(1,var)*E_ins/Cond_h(1,var)+...
                     +(1/(E_cbl*SC_w(1,var)/SC_h(1,var))+2/(E_jckt*Cond_w(1,var)/JT(1,var))+2/(E_ins*Cond_w(1,var)/INS(1,var)))^-1;
Ke_cavo_tor(1,var) = 2*E_jckt*JT(1,var)/Cond_w(1,var)+2*INS(1,var)*E_ins/Cond_w(1,var)+...
                     +(1/(E_cbl*SC_h(1,var)/SC_w(1,var))+2/(E_jckt*JT(1,var)/Cond_w(1,var))+2/(E_ins*INS(1,var)/Cond_w(1,var)))^-1;
K_jckt = 2*JT(1,var)/Cond_h(1,var)*E_jckt;
dcr_jckt = K_jckt/Ke_cavo_rad(1,var);
r_steel = (Cond_w(n_layers)-2*INS(n_layers))/(2*JT(n_layers));
S_rm = p_rs*r_steel*scf*dcr_jckt; % Radial memebrane stress Jacket innermost layer  

% if cell2sym(type_cable(1,var)) == 'HTS'
%    if S_rm_Core = (p_rs/2+S_z_JT/2)> S_amm_Cu
%        continue
%    end
% end
%% Cechk membrane stress on Jacket 
clear Rk_ DTF S_T_JT S_T_VT Ke_WP_rad Ke_WP_tor K_ps_rad K_ps_tor dcr_vault_tor dcr_vault_rad...
                    K_vault_rad K_vault_tor K_lat_rad K_lat_tor
DTF = 0.05; % Vaul width 
S_T_JT = 1e30;
S_T_VT = 1e30; % Test value for vault Tresca stress
while S_T_VT > S_amm_VT || S_c_VT > S_amm_VT/1.3 || S_T_JT > S_amm_JT || S_rm_JT > S_amm_JT/1.3
    DTF = DTF+0.001; % Se supero Tresca incremento spessore naso TF
    Rk_ = Rj_-DTF; % Innermost Case radius 
    CASE_w_l = 2*Rk_*tan(theta_TF/2); % Case low part width 
    A_tot = (2*Ri_*tan(theta_TF/2)+2*Rk_*tan(theta_TF/2))*(Ri_-Rk_)/2;
    A_CASE = A_tot-A_WP;
    A_VT = (2*Rj_*tan(theta_TF/2)+2*Rk_*tan(theta_TF/2))*(Rj_-Rk_)/2; % Vault section trp
    A_VT_circ = pi*(Rj_^2-Rk_^2)*1/n_TF; % Vault section circ
    % Rigidezze    
    Ke_WP_rad = sum(1./(Ke_cavo_rad(1:n_layers).*n_turns(1:n_layers)))^-1;    
    Ke_WP_tor = sum((n_turns(1:n_layers)./Ke_cavo_tor(1:n_layers)).^-1);   
    
    K_ps_rad = E_case/dr_plasma_side*(2*Ri_*tan(theta_TF/2));
    K_ps_tor = (E_case*dr_plasma_side/(2*Ri_*tan(theta_TF/2)));
    K_lat_rad = 2*E_case*lateral_w/WP_h;
    K_lat_tor = (2/(E_case*WP_h/lateral_w))^-1;
    K_vault_rad = (E_case*Rk_*2*tan(theta_TF/2)/DTF);
    K_vault_tor = (E_case*DTF/(Rk_*2*tan(theta_TF/2)));

%     Ke_case_rad = (1/(Ke_WP_rad+2*K_lat_rad)+1/K_vault_rad)^-1;
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
    S_c_VT = 2/(1-beta^2)*p_rs*dcr_vault_tor*1.1*A_VT/A_VT_circ; % @ Rk_
    S_l_VT = -1/(1-beta^2)*p_rs; % @ Rk_    

    k_bf = 0.5*log(RTFo(dp)/RTFi(dp)); % k bending free
    T_bf = 0.5*(k_bf*n_TF*(n_spire_(1,1)*Iop)^2*Mu_0/(2*pi)); % Hoop tension along TF longitudinal axis

    S_z =  T_bf/(A_J_tot+A_CASE)*1.1; % Normal tension on Case stell section 
    S_T_VT = (S_z+S_c_VT)*1.1; % Vault Tresca stress 
    S_T_JT = (S_z+S_rm_JT);  % Jacket Tresca stress      
end  
V = L*Iop/Tau_discharge*1e-3; % Tensione max singolo TF
%%
% figure('units','normalized','outerposition',[0.25 0.1 0.5 0.9])
figure('units','normalized','outerposition',[0 0 1 1])
% c = uisetcolor([1 1 0],'Select a color')
map = jet(n_layers); 
J_SC = Iop./EQV_Cable_Area(1:n_layers).*1e-6;
fill([-CASE_w_l/2 -CASE_w/2 CASE_w/2 CASE_w_l/2 -CASE_w_l/2],[Rk_ Ri_ Ri_ Rk_ Rk_],[0.7 0.7 0.7]), hold on
for var = 1:n_layers
    for t = 1:n_turns(var)                                
        % map(var,:) = [min(J_SC)/J_SC(var)    min(J_SC)/J_SC(var)^2    min(J_SC)/J_SC(var)];  
        rectangle('Position',[(-Cond_w(1,var)-Cond_w(1,var)*n_turns(var)/2)+Cond_w(1,var)*t Ri(1,var) Cond_w(1,var) Cond_h(1,var)],...
            'FaceColor',[0.0588    1.0000    1.0000])
        hold on
        rectangle('Position',[(-Cond_w(1,var)-Cond_w(1,var)*n_turns(var)/2+INS(1,var))+Cond_w(1,var)*t Ri(1,var)+INS(1,var) Cond_w(1,var)-2*INS(1,var) Cond_h(1,var)-2*INS(1,var)],...
            'FaceColor',[0.5    0.5    0.5])
        hold on                     
        if shape_cable == 'RIS_'                      
            rectangle('Position',[(Cond_w(1,var)*(-1-n_turns(var)/2)+Cond_w(1,var)/2-SC_w(1,var)/2)+Cond_w(1,var)*t Ri(1,var)+JT(1,var)+2*INS(1,var)-2*INS_grades SC_w(1,var) SC_w(1,var)],...
            'FaceColor',map(n_layers-var+1,:),'Curvature',[1 1])          
        else
            rectangle('Position',[(-Cond_w(1,var)-Cond_w(1,var)*n_turns(var)/2+JT(1,var)+INS(1,var))+Cond_w(1,var)*t Ri(1,var)+JT(1,var)+INS(1,var) Cond_w(1,var)-2*JT(1,var)-2*INS(1,var) Cond_h(1,var)-2*JT(1,var)-2*INS(1,var)],...
            'FaceColor',map(n_layers-var+1,:),'Curvature',[0 0])
        end         
    end
end
axis equal
set(gca,'FontSize',20);
colormap(map);
colorbar
c = colorbar;
c.Label.FontSize = 20;
c.Label.Color = 'k';
c.Label.Rotation = 90;
c.Label.String = 'B $[T]$';
c.Label.Interpreter = 'latex';
c.TickLabelInterpreter = 'latex';
mod_max = round(max(B_layers),1);
mod_min = round(min(B_layers),1);
clim([mod_min  mod_max]);
MAXMIN=get(c,'Limits');
T = round(linspace(MAXMIN(1),MAXMIN(2),8),2);
set(c,'Ticks',T)
% TL = arrayfun(@(x) sprintf('%.2f',x),T,'un',0);
xlabel('[m]','fontsize',20,'Interpreter','latex')
ylabel('[m]','fontsize',20,'Interpreter','latex')

TitleName2 = sprintf('n grades=%G, N=%G , Iop=%G kA, Sa(VT)=%G MPa, Sa(JT)=%G MPa, tau=%G s',n_grades,n_spire_(1), round(Iop*1e-3),round(S_amm_VT*1e-6),round(S_amm_JT*1e-6),round(Tau_discharge));
title(TitleName2,'fontsize',20,'Interpreter','latex');

ylim([Rk_*0.8 Ri_*1.2])
xlim([-CASE_w*2 CASE_w*2])
dx_ = -CASE_w*1.01;
plot([dx_ dx_],[0 Rk_],'--.r','linewidth',1); hold on;  
plot([dx_ dx_],[Rk_ Rj_],'--.g','linewidth',1); hold on; 
plot([dx_ dx_],[Rj_ Ri_],'--.b','linewidth',1); hold on;  
plot([-dx_ -dx_],[Rj_ (Rj_+WP_h)],'--.k','linewidth',1); hold on;  
plot([-dx_ -dx_],[Rk_ Rj_],'--.k','linewidth',1); hold on; 
plot([-dx_ -dx_],[(Rj_+WP_h) Ri_],'--.k','linewidth',1); hold on; 

plot([-WP_w/2 WP_w/2],[Ri_*1.03 Ri_*1.03],'--.k','linewidth',1); hold on;  
plot([-CASE_w/2 CASE_w/2],[Ri_*1.08 Ri_*1.08],'--.k','linewidth',1); hold on; 

dx_ = -CASE_w;
text(dx_,Rk_,sprintf('Rk=%G m',round(Rk_,2)),'fontsize',12,'Interpreter','latex')
text(dx_,Rj_,sprintf('Rj=%G m',round(Rj_,2)),'fontsize',12,'Interpreter','latex')
text(dx_,Ri_,sprintf('Ri=%G m',round(Ri_,2)),'fontsize',12,'Interpreter','latex')
text(-dx_*1.03,Rj_+WP_h/2,sprintf('WP(h)=%G m',round(WP_h,2)),'fontsize',12,'Interpreter','latex')
text(-dx_*1.03,(Rk_+Rj_)/2,sprintf('Nose(h)=%G m',round(Rj_-Rk_,2)),'fontsize',12,'Interpreter','latex')
text(-dx_*1.03,((Rj_+WP_h)+Ri_)/2,sprintf('dps(h)=%G m',round(dr_plasma_side,2)),'fontsize',12,'Interpreter','latex')
text(0,Ri_*1.1,sprintf('Case(w)=%G m',round(CASE_w,2)),'HorizontalAlignment','center','fontsize',12,'Interpreter','latex')
text(0,Ri_*1.05,sprintf('WP(w)=%G m',round(WP_w,2)),'HorizontalAlignment','center','fontsize',12,'Interpreter','latex')

if printcond == 1, print('-dpng','-r300',TitleName2), end
%% Plot WP
% c = uisetcolor([1 1 0],'Select a color')
% figure('units','normalized','outerposition',[0.25 0.1 0.5 0.9])
figure('units','normalized','outerposition',[0 0 1 1])
J_eng = Iop./(Cond_w.*Cond_h).*1e-6;
J_CICC = Iop./(EQV_Cable_Area(1:n_layers)).*1e-6;
mod_max = max(J_eng(1:n_layers));
mod_min = min(J_eng(1:n_layers));
if mod_max == mod_min
    clim([0 mod_max]);
    map = jet(1); 
    colorvar = ones(1,n_layers);
else
    clim([mod_min mod_max]);
    map = jet(n_grades); 
    colorvar = linspace(1,n_grades,n_grades);
end
fill([-CASE_w_l/2 -CASE_w/2 CASE_w/2 CASE_w_l/2 -CASE_w_l/2],[Rk_ Ri_ Ri_ Rk_ Rk_],[0.7 0.7 0.7]), hold on
n_layerspergrade = grading_at; n_layerspergrade(end) = n_layers; 
for varg = 1:n_grades
    for var = n_layerspergrade(varg):n_layerspergrade(varg+1)  
        for t = 1:n_turns(var)                           
            % map(var,:) = [min(J_SC)/J_SC(var)    min(J_SC)/J_SC(var)^2    min(J_SC)/J_SC(var)];  
            rectangle('Position',[(Cond_w(1,var)*(-1-n_turns(var)/2))+Cond_w(1,var)*t Ri(1,var) Cond_w(1,var) Cond_h(1,var)],...
                'FaceColor',map(colorvar(varg),:))
            hold on
            rectangle('Position',[(Cond_w(1,var)*(-1-n_turns(var)/2)+INS(1,var))+Cond_w(1,var)*t Ri(1,var)+INS(1,var) Cond_w(1,var)-2*INS(1,var) Cond_h(1,var)-2*INS(1,var)],...
                'FaceColor',map(colorvar(varg),:))
           hold on
            if shape_cable == 'RIS_'                      
                rectangle('Position',[(Cond_w(1,var)*(-1-n_turns(var)/2)+Cond_w(1,var)/2-SC_w(1,var)/2)+Cond_w(1,var)*t Ri(1,var)+JT(1,var)+2*INS(1,var)-2*INS_grades SC_w(1,var) SC_w(1,var)],...
                'FaceColor',map(colorvar(varg),:),'Curvature',[1 1])          
            else
                rectangle('Position',[(-Cond_w(1,var)-Cond_w(1,var)*n_turns(var)/2+JT(1,var)+INS(1,var))+Cond_w(1,var)*t Ri(1,var)+JT(1,var)+INS(1,var) Cond_w(1,var)-2*JT(1,var)-2*INS(1,var) Cond_h(1,var)-2*JT(1,var)-2*INS(1,var)],...
                'FaceColor',map(colorvar(varg),:),'Curvature',[0 0])
            end          
        end
    end
end
set(gca,'FontSize',20);
axis equal
colormap(map);
colorbar
c = colorbar;
c.Label.FontSize = 20;
c.Label.Color = 'k';
c.Label.Rotation = 90;
c.Label.String = 'J eng $[A/mm^2]$';
c.Label.Interpreter = 'latex';
c.TickLabelInterpreter = 'latex';
mod_max = max(J_eng(1:n_layers));
mod_min = min(J_eng(1:n_layers));
if mod_max == mod_min
    clim([0 mod_max]);
    map = jet(1); 
    colorvar = ones(1,n_layers);
    divis = 1;
else
    clim([mod_min mod_max]);
    map = jet(n_grades); 
    divis = n_grades+1;
end
MAXMIN=get(c,'Limits');
T = floor(linspace(MAXMIN(1),MAXMIN(2),divis));
set(c,'Ticks',T)
% TL = arrayfun(@(x) sprintf('%.2f',x),T,'un',0);
xlabel('[m]','fontsize',20,'Interpreter','latex')
ylabel('[m]','fontsize',20,'Interpreter','latex')
TitleName1 = sprintf('n grades=%G, N=%G , Iop=%G kA, Sa(VT)=%G MPa, Sa(JT)=%G MPa, tau=%G s ',n_grades,n_spire_(1), round(Iop*1e-3),round(S_amm_VT*1e-6),round(S_amm_JT*1e-6),round(Tau_discharge));
title(TitleName1,'fontsize',20,'Interpreter','latex');
ylim([Rk_*0.8 Ri_*1.2])
xlim([-CASE_w*2 CASE_w*2])
dx_ = -CASE_w*1.01;
plot([dx_ dx_],[0 Rk_],'--.r','linewidth',1); hold on;  
plot([dx_ dx_],[Rk_ Rj_],'--.g','linewidth',1); hold on; 
plot([dx_ dx_],[Rj_ Ri_],'--.b','linewidth',1); hold on;  
plot([-dx_ -dx_],[Rj_ (Rj_+WP_h)],'--.k','linewidth',1); hold on;  
plot([-dx_ -dx_],[Rk_ Rj_],'--.k','linewidth',1); hold on; 
plot([-dx_ -dx_],[(Rj_+WP_h) Ri_],'--.k','linewidth',1); hold on; 

plot([-WP_w/2 WP_w/2],[Ri_*1.03 Ri_*1.03],'--.k','linewidth',1); hold on;  
plot([-CASE_w/2 CASE_w/2],[Ri_*1.08 Ri_*1.08],'--.k','linewidth',1); hold on; 

dx_ = -CASE_w;
text(dx_,Rk_,sprintf('Rk=%G m',round(Rk_,2)),'fontsize',12,'Interpreter','latex')
text(dx_,Rj_,sprintf('Rj=%G m',round(Rj_,2)),'fontsize',12,'Interpreter','latex')
text(dx_,Ri_,sprintf('Ri=%G m',round(Ri_,2)),'fontsize',12,'Interpreter','latex')
text(-dx_*1.03,Rj_+WP_h/2,sprintf('WP(h)=%G m',round(WP_h,2)),'fontsize',12,'Interpreter','latex')
text(-dx_*1.03,(Rk_+Rj_)/2,sprintf('Nose(h)=%G m',round(Rj_-Rk_,2)),'fontsize',12,'Interpreter','latex')
text(-dx_*1.03,((Rj_+WP_h)+Ri_)/2,sprintf('dps(h)=%G m',round(dr_plasma_side,2)),'fontsize',12,'Interpreter','latex')
text(0,Ri_*1.1,sprintf('Case(w)=%G m',round(CASE_w,2)),'HorizontalAlignment','center','fontsize',12,'Interpreter','latex')
text(0,Ri_*1.05,sprintf('WP(w)=%G m',round(WP_w,2)),'HorizontalAlignment','center','fontsize',12,'Interpreter','latex')

if printcond == 1, print('-dpng','-r300',TitleName1), end
%% Plot WP - Jeng comparison
% 
% figure('units','normalized','outerposition',[0 0 1 1])
% jeng_ = ceil(Iop./(Cond_w(1).*Cond_h(1)).*1e-6);
% J_CICC = Iop./(EQV_Cable_Area(1:n_layers)).*1e-6;
% mod_max = 26;
% mod_min = 8;
% map = jet(mod_max-mod_min); 
% colorvar = linspace(1,mod_max-mod_min,mod_max-mod_min);
% 
% fill([-CASE_w_l/2 -CASE_w/2 CASE_w/2 CASE_w_l/2 -CASE_w_l/2],[Rk_ Ri_ Ri_ Rk_ Rk_],[0.7 0.7 0.7]), hold on
% for var = 1:n_layers  
%     for t = 1:n_turns(var)                           
%         % map(var,:) = [min(J_SC)/J_SC(var)    min(J_SC)/J_SC(var)^2    min(J_SC)/J_SC(var)];  
%         rectangle('Position',[(Cond_w(1,var)*(-1-n_turns(var)/2))+Cond_w(1,var)*t Ri(1,var) Cond_w(1,var) Cond_h(1,var)],...
%             'FaceColor',map(jeng_-mod_min,:))
%         hold on
%         rectangle('Position',[(Cond_w(1,var)*(-1-n_turns(var)/2)+INS(1,var))+Cond_w(1,var)*t Ri(1,var)+INS(1,var) Cond_w(1,var)-2*INS(1,var) Cond_h(1,var)-2*INS(1,var)],...
%             'FaceColor',map(jeng_-mod_min,:))
%        hold on
%         if shape_cable == 'RIS_'                      
%             rectangle('Position',[(Cond_w(1,var)*(-1-n_turns(var)/2)+Cond_w(1,var)/2-SC_w(1,var)/2)+Cond_w(1,var)*t Ri(1,var)+JT(1,var)+2*INS(1,var)-2*INS_grades SC_w(1,var) SC_w(1,var)],...
%             'FaceColor',map(jeng_-mod_min,:),'Curvature',[1 1])          
%         else
%             rectangle('Position',[(-Cond_w(1,var)-Cond_w(1,var)*n_turns(var)/2+JT(1,var)+INS(1,var))+Cond_w(1,var)*t Ri(1,var)+JT(1,var)+INS(1,var) Cond_w(1,var)-2*JT(1,var)-2*INS(1,var) Cond_h(1,var)-2*JT(1,var)-2*INS(1,var)],...
%             'FaceColor',map(jeng_-mod_min,:),'Curvature',[0 0])
%         end          
%     end
% end
% set(gca,'FontSize',20);
% axis equal
% colormap(map);
% colorbar
% c = colorbar;
% c.Label.FontSize = 20;
% c.Label.Color = 'k';
% c.Label.Rotation = 90;
% c.Label.String = 'J eng $[A/mm^2]$';
% c.Label.Interpreter = 'latex';
% c.TickLabelInterpreter = 'latex';
% clim([mod_min mod_max]);
% MAXMIN=get(c,'Limits');
% T = floor(linspace(MAXMIN(1),MAXMIN(2),9));
% set(c,'Ticks',T)
% % TL = arrayfun(@(x) sprintf('%.2f',x),T,'un',0);
% xlabel('[m]','fontsize',20,'Interpreter','latex')
% ylabel('[m]','fontsize',20,'Interpreter','latex')
% TitleName = sprintf('N=%G , Iop=%G kA, Sy(VT)=%G MPa, Sy(JT)=%G MPa, tau=%G s ',n_spire_(1), round(Iop*1e-3),round(S_amm_VT*1e-6),round(S_amm_JT*1e-6),round(Tau_discharge));
% title(TitleName,'fontsize',20,'Interpreter','latex');
% ylim([Rk_*0.8 Ri_*1.2])
% xlim([-CASE_w*2 CASE_w*2])
% dx_ = -CASE_w*1.01;
% plot([dx_ dx_],[0 Rk_],'--.r','linewidth',1); hold on;  
% plot([dx_ dx_],[Rk_ Rj_],'--.g','linewidth',1); hold on; 
% plot([dx_ dx_],[Rj_ Ri_],'--.b','linewidth',1); hold on;  
% plot([-dx_ -dx_],[Rj_ (Rj_+WP_h)],'--.k','linewidth',1); hold on;  
% plot([-dx_ -dx_],[Rk_ Rj_],'--.k','linewidth',1); hold on; 
% plot([-dx_ -dx_],[(Rj_+WP_h) Ri_],'--.k','linewidth',1); hold on; 
% 
% plot([-WP_w/2 WP_w/2],[Ri_*1.03 Ri_*1.03],'--.k','linewidth',1); hold on;  
% plot([-CASE_w/2 CASE_w/2],[Ri_*1.08 Ri_*1.08],'--.k','linewidth',1); hold on; 
% 
% dx_ = -CASE_w;
% text(dx_,Rk_,sprintf('Rk=%G m',round(Rk_,2)),'fontsize',12,'Interpreter','latex')
% text(dx_,Rj_,sprintf('Rj=%G m',round(Rj_,2)),'fontsize',12,'Interpreter','latex')
% text(dx_,Ri_,sprintf('Ri=%G m',round(Ri_,2)),'fontsize',12,'Interpreter','latex')
% text(-dx_*1.03,Rj_+WP_h/2,sprintf('WP(h)=%G m',round(WP_h,2)),'fontsize',12,'Interpreter','latex')
% text(-dx_*1.03,(Rk_+Rj_)/2,sprintf('Nose(h)=%G m',round(Rj_-Rk_,2)),'fontsize',12,'Interpreter','latex')
% text(-dx_*1.03,((Rj_+WP_h)+Ri_)/2,sprintf('dps(h)=%G m',round(dr_plasma_side,2)),'fontsize',12,'Interpreter','latex')
% text(0,Ri_*1.1,sprintf('Case(w)=%G m',round(CASE_w,2)),'HorizontalAlignment','center','fontsize',12,'Interpreter','latex')
% text(0,Ri_*1.05,sprintf('WP(w)=%G m',round(WP_w,2)),'HorizontalAlignment','center','fontsize',12,'Interpreter','latex')
% text(0,Ri_*1.05,sprintf('WP(w)=%G m',round(WP_w,2)),'HorizontalAlignment','center','fontsize',12,'Interpreter','latex')
% 
% if printcond == 1, print('-dpng','-r300',TitleName1), end
%% Plot WP - B comparison
% 
% figure('units','normalized','outerposition',[0 0 1 1])
% jeng_ = floor(B_PHI_TF);
% mod_max = 12.5;
% mod_min = 0;
% map = jet(ceil(mod_max-mod_min));
% colorvar = linspace(1,mod_max-mod_min,mod_max-mod_min);
% 
% fill([-CASE_w_l/2 -CASE_w/2 CASE_w/2 CASE_w_l/2 -CASE_w_l/2],[Rk_ Ri_ Ri_ Rk_ Rk_],[0.7 0.7 0.7]), hold on
% for var = 1:n_layers  
%     for t = 1:n_turns(var)                           
%         % map(var,:) = [min(J_SC)/J_SC(var)    min(J_SC)/J_SC(var)^2    min(J_SC)/J_SC(var)];  
%         rectangle('Position',[(Cond_w(1,var)*(-1-n_turns(var)/2))+Cond_w(1,var)*t Ri(1,var) Cond_w(1,var) Cond_h(1,var)],...
%             'FaceColor',[0.5    0.5    0.5])
%         hold on
%         rectangle('Position',[(Cond_w(1,var)*(-1-n_turns(var)/2)+INS(1,var))+Cond_w(1,var)*t Ri(1,var)+INS(1,var) Cond_w(1,var)-2*INS(1,var) Cond_h(1,var)-2*INS(1,var)],...
%             'FaceColor',[0.5    0.5    0.5])
%        hold on
%         if shape_cable == 'RIS_'                      
%             rectangle('Position',[(Cond_w(1,var)*(-1-n_turns(var)/2)+Cond_w(1,var)/2-SC_w(1,var)/2)+Cond_w(1,var)*t Ri(1,var)+JT(1,var)+2*INS(1,var)-2*INS_grades SC_w(1,var) SC_w(1,var)],...
%             'FaceColor',map(jeng_-mod_min,:),'Curvature',[1 1])          
%         else
%             rectangle('Position',[(-Cond_w(1,var)-Cond_w(1,var)*n_turns(var)/2+JT(1,var)+INS(1,var))+Cond_w(1,var)*t Ri(1,var)+JT(1,var)+INS(1,var) Cond_w(1,var)-2*JT(1,var)-2*INS(1,var) Cond_h(1,var)-2*JT(1,var)-2*INS(1,var)],...
%             'FaceColor',map(jeng_-mod_min-var+1,:),'Curvature',[0 0])
%         end          
%     end
% end
% set(gca,'FontSize',20);
% axis equal
% colormap(map);
% colorbar
% c = colorbar;
% c.Label.FontSize = 20;
% c.Label.Color = 'k';
% c.Label.Rotation = 90;
% c.Label.String = 'B $[T]$';
% c.Label.Interpreter = 'latex';
% c.TickLabelInterpreter = 'latex';
% clim([mod_min mod_max]);
% MAXMIN=get(c,'Limits');
% T = floor(linspace(MAXMIN(1),MAXMIN(2),9));
% set(c,'Ticks',T)
% % TL = arrayfun(@(x) sprintf('%.2f',x),T,'un',0);
% xlabel('[m]','fontsize',20,'Interpreter','latex')
% ylabel('[m]','fontsize',20,'Interpreter','latex')
% TitleName = sprintf('N=%G , Iop=%G kA, Sy(VT)=%G MPa, Sy(JT)=%G MPa, tau=%G s ',n_spire_(1), round(Iop*1e-3),round(S_amm_VT*1e-6),round(S_amm_JT*1e-6),round(Tau_discharge));
% title(TitleName,'fontsize',20,'Interpreter','latex');
% ylim([Rk_*0.8 Ri_*1.2])
% xlim([-CASE_w*2 CASE_w*2])
% dx_ = -CASE_w*1.01;
% plot([dx_ dx_],[0 Rk_],'--.r','linewidth',1); hold on;  
% plot([dx_ dx_],[Rk_ Rj_],'--.g','linewidth',1); hold on; 
% plot([dx_ dx_],[Rj_ Ri_],'--.b','linewidth',1); hold on;  
% plot([-dx_ -dx_],[Rj_ (Rj_+WP_h)],'--.k','linewidth',1); hold on;  
% plot([-dx_ -dx_],[Rk_ Rj_],'--.k','linewidth',1); hold on; 
% plot([-dx_ -dx_],[(Rj_+WP_h) Ri_],'--.k','linewidth',1); hold on; 
% 
% plot([-WP_w/2 WP_w/2],[Ri_*1.03 Ri_*1.03],'--.k','linewidth',1); hold on;  
% plot([-CASE_w/2 CASE_w/2],[Ri_*1.08 Ri_*1.08],'--.k','linewidth',1); hold on; 
% 
% dx_ = -CASE_w;
% text(dx_,Rk_,sprintf('Rk=%G m',round(Rk_,2)),'fontsize',12,'Interpreter','latex')
% text(dx_,Rj_,sprintf('Rj=%G m',round(Rj_,2)),'fontsize',12,'Interpreter','latex')
% text(dx_,Ri_,sprintf('Ri=%G m',round(Ri_,2)),'fontsize',12,'Interpreter','latex')
% text(-dx_*1.03,Rj_+WP_h/2,sprintf('WP(h)=%G m',round(WP_h,2)),'fontsize',12,'Interpreter','latex')
% text(-dx_*1.03,(Rk_+Rj_)/2,sprintf('Nose(h)=%G m',round(Rj_-Rk_,2)),'fontsize',12,'Interpreter','latex')
% text(-dx_*1.03,((Rj_+WP_h)+Ri_)/2,sprintf('dps(h)=%G m',round(dr_plasma_side,2)),'fontsize',12,'Interpreter','latex')
% text(0,Ri_*1.1,sprintf('Case(w)=%G m',round(CASE_w,2)),'HorizontalAlignment','center','fontsize',12,'Interpreter','latex')
% text(0,Ri_*1.05,sprintf('WP(w)=%G m',round(WP_w,2)),'HorizontalAlignment','center','fontsize',12,'Interpreter','latex')
% text(0,Ri_*1.05,sprintf('WP(w)=%G m',round(WP_w,2)),'HorizontalAlignment','center','fontsize',12,'Interpreter','latex')
% if printcond == 1, print('-dpng','-r300',TitleName), end
%% FEM

% ! "C:\Program Files\ANSYS Inc\v221\ansys\bin\winx64\ANSYS221.exe" -b -np 6 -i MAG.dat -o MAG.out -j MAG
% 
% delete('mag.out'); delete('str.out')
% for i=0:13
%     delete(sprintf('MAG%G.err',i));     delete(sprintf('STR%G.err',i));
%     delete(sprintf('MAG%G.out',i));     delete(sprintf('STR%G.out',i));
% end


