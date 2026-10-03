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
n_turns = ones(1,n_layers)*[20];
% n_turns = [10	10	10	10	10	10	10	10	10	10	10	10	10	10	10];
n_grades = 1;
grading_at = [1 (n_layers+1)]; %% il grade viene applicato a partire dal layer indicato (incluso)
%%
lateral_w = 0.20;
min_JT = 0.002; 
V_MAX = 10e3; % [V]
%%
dp = 1; 
%% input da esplorazione 
n_TF = 16;  
B0(dp) = 4.389;       
R0(dp) = 8.6;  
A(dp) = 2.8;
wb = 1.82;
RTFi(dp) = R0(dp)-R0(dp)/A(dp)-wb;    % [m]
ripple = 0.01; 
RTFo(dp) = (R0(dp)+R0(dp)/A(dp))*(1/ripple)^(1/n_TF)+1;    % [m]
R_VV = RTFi(dp)*1.05;                  % [m]
%%
corr_B_WP   = 1.08;
WP_SC_type  = 100;   % full_LTS = 100 - full_HTS = 101 - Hybrid = 102
shape_cable = 200;   % 200 = RIS - 201 = Rect
%% Ammissibili acciaio
S_amm_JT                   = (667*1e6);
S_amm_VT                   = (667*1e6);
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
E_ins = 10;
E_cbl_LTS = 0;
E_case = 205;
%% WP dimensioning
CASE_w          = 2*R_TF_Innerleg*tan(theta_TF/2); % Case width
dr_plasma_side  = 0.06; % [m] Spessore case fronte plasma
GoundIns        = 0.008; % Ground insulation WP 
INS_grades      = 0.0005; % Spessore isolante tra i Grades
%% Dati configurazioni da iterare 
iter = iter+1;                                       
%% Turns in each grade
n_spire_ = zeros(1,n_layers);
n_spire_(1,1) = sum(n_turns);
for var = 2:n_layers 
    n_spire_(1,var) = n_spire_(1,var-1)-n_turns(1,var-1);
end   
%% Preallocations
Cond_h = zeros(1,maxdim); 
Cond_w = zeros(1,maxdim); 
clear S_Cable;
% S_Cable = zeros(1,maxdim); 
JT = zeros(1,maxdim); 
r_cable = zeros(1,maxdim); 
type_cable = cell(1, maxdim);
for k = 1:maxdim
    type_cable(1,k) = {'---'};
end 
S_Cable = zeros(1,n_layers); N_tot = zeros(1,maxdim);
N_Sc = zeros(1,maxdim); Ic_sc = zeros(1,n_layers);
N_sc = zeros(1,maxdim); N_Cu = zeros(1,maxdim);
THS = zeros(1,maxdim); N_sc0 = zeros(1,n_layers);
r_SC = zeros(1,n_layers); S_Cu_HTS = zeros(1,maxdim);
S_REBCO = zeros(1,maxdim); SC_w = zeros(1,n_layers);
R_J = zeros(1,n_layers); SC_h = zeros(1,n_layers);
S_CICC = zeros(1,n_layers); S_JT = zeros(1,n_layers);
Ri = zeros(1,n_layers); Re = zeros(1,n_layers);        
%% Operative current definition
Iop = ceil(NI/n_spire_(1,1));   
%% Induttanza TF modello shell
clear L Tau_discharge
% L = Mu_0*R0(dp)*(n_TF*n_spire_(1,1))^2*(1-sqrt(1-(R0(dp)-RTFi(dp))/R0(dp)))/n_TF*1.1;
% L = Mu_0*R0(dp)*(n_TF*n_spire_(1,1))^2*(1-sqrt(1-(RTFo(dp)-RTFi(dp))/2/R0(dp)))/n_TF;

R_WP_IL = R_TF_Innerleg;
R_WP_OL = R_TF_Outerleg;
k_bf = 0.5*log(R_WP_OL/R_WP_IL); % k bending free
r_bf = sqrt(RTFo(dp)*RTFi(dp));
Ntot = (n_TF*n_spire_(1,1));
L = Mu_0*(Ntot*k_bf)^2*r_bf/2*(besseli(0, k_bf)+ 2*besseli(1, k_bf)+ besseli(2, k_bf))/n_TF;

Tau_discharge2 = (L*Iop/V_MAX); % [s] - scarico in gruppi di tre bobine
Tau_discharge = max([Tau_discharge1 Tau_discharge2 4]);
E = 1/2*L*n_TF*Iop^2*1e-6;

%% Definition of the magnetic field peak in each grade (linera beahvior from B(Re)= Bmax to B(Ri) = 0                 
Re(1,1) = R_TF_Innerleg-dr_plasma_side-GoundIns; % WP innerl-leg outer radius (non-insulated)
B_TF = B_PHI_TF; % (n_TF*n_spire_(1,1)*Iop*Mu_0)/(2*pi*Re(1,1)); 
B_layers = B_TF.*(n_spire_./n_spire_(1,1)); % Ottengo B considerando un andamento lineare nel WP
%%
for g = 1:n_grades
    var =  grading_at(g);
    B = B_layers(1,var);       
    [type_cable(1,var),N_Cu(1,var),N_Sc(1,var),N_tot(1,var),S_Cable(1,var),S_REBCO(1,var),S_Cu_HTS(1,var),THS(1,var)] ...
    = CICC_DEMO(B,Iop,Tau_discharge,WP_SC_type); 

    type_cable(1,var:n_layers) = type_cable(1,var);
    N_Cu(1,var:n_layers) = N_Cu(1,var);
    N_Sc(1,var:n_layers) = N_Sc(1,var);
    N_tot(1,var:n_layers) = N_tot(1,var);
    S_Cable(1,var:n_layers) = S_Cable(1,var);
    S_REBCO(1,var:n_layers) = S_REBCO(1,var);
    S_Cu_HTS(1,var:n_layers) = S_Cu_HTS(1,var);
    THS(1,var:n_layers) = THS(1,var);
end

%%  Definisco sezioni cavi in ogni layers 
R_WP_IL = R_TF_Innerleg;
R_WP_OL = R_TF_Outerleg;
k_bf = 0.5*log(R_WP_OL/R_WP_IL); % k bending free
T_bf = 0.5*(k_bf*n_TF*(n_spire_(1,1)*Iop)^2*Mu_0/(2*pi)); % Hoop tension along TF longitudinal axis
T_bf_  = T_bf*1e-6;
% Fr = (n_TF*(n_spire_(1,1)*Iop)^2*Mu_0/2)*(1-(1/sqrt(1-(RTFi(dp)/R0(dp))^2)));
F0 = n_TF*(n_spire_(1,1)*Iop)^2*Mu_0/4/pi;
H = 1;
Fr = -F0*H/RTFi(dp)*1e-6;

WP_w0(1,1) = 2*Re(1,1)*tan(theta_TF/2)-lateral_w*2-GoundIns*2; % massimo ingombro toroidale WP
S_z_JT =  T_bf/(WP_w0(1,1)^2)/2;
p_rs = B_TF^2/(2*Mu_0);                                 % Magnetic pressure, thin WP

for var = 1:n_layers 
    % Definisco dimensioni interne del cavo CICC per l'i-esimo grade
    r_SC(1,var)              = (0.005)*Increm; %JT(1,var);   % Impogno raggio curvatura corner cavo = a spessore Jakcet
    tins(1,var)              = (0.001)*Increm; % Isolante di spira
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
        JT(1,var) = min_JT;
        S_rm_grades(1,var) = 1e30;
        while S_rm_grades(1,var) > S_amm_JT/1.3
            JT(1,var)                = JT(1,var) + 1e-4;
            SC_w(1,var)              = Cond_w(1,var)-JT(1,var)*2-tins(1,var)*2; % Ottengo larghezza cavo SC 
            SC_h(1,var)              = (S_Cable(1,var)+(4-pi)*r_SC(1,var)^2)/SC_w(1,var); % Ottengo altezza cavo SC
            R_J(1,var)               = r_SC(1,var) + JT(1,var); % Raggio di curvatura corner Jacket
            Cond_h(1,var)            = SC_h(1,var) + JT(1,var)*2+tins(1,var)*2; % Ottengo l'altezza del cavo Rect                        
            Ke_cavo_rad(1,var)       = 2*E_jckt*JT(1,var)/Cond_h(1,var)+2*tins(1,var)*E_ins/Cond_h(1,var)+...
                                        +(1/(E_cbl*SC_w(1,var)/SC_h(1,var))+2/(E_jckt*Cond_w(1,var)/JT(1,var))+2/(E_ins*Cond_w(1,var)/tins(1,var)))^-1;
            Ke_cavo_tor(1,var)       = 2*E_jckt*JT(1,var)/Cond_w(1,var)+2*tins(1,var)*E_ins/Cond_w(1,var)+...
                                        +(1/(E_cbl*SC_h(1,var)/SC_w(1,var))+2/(E_jckt*JT(1,var)/Cond_w(1,var))+2/(E_ins*tins(1,var)/Cond_w(1,var)))^-1;
            K_jckt                   = 2*JT(1,var)/Cond_h(1,var)*E_jckt;
            dcr_jckt                 = K_jckt/Ke_cavo_rad(1,var);
            r_steel                  = (Cond_w(1,var))/(2*JT(1,var));
            S_rm_grades(1,var)       = p_rs*r_steel*dcr_jckt+S_z_JT;      % Radial memebrane stress Jacket innermost layer 
        end  
    end
    S_CICC(1,var)            = ((Cond_w(1,var)-2*tins(1,var))*(Cond_h(1,var)-2*tins(1,var)))-((4-pi)*R_J(1,var)^2); % Sezione cavo non isolato
    S_JT(1,var)              = S_CICC(1,var)-S_Cable(1,var); % Sezione di acciaio                                
    Ri(1,var)                = Re(1,var)-Cond_h(1,var); % Raggio interno grade i-esimo del WP
    Re(1,var+1)              = Ri(1,var)-INS_grades;   % Raggio esterno grade (i+1)-esimo del WP     
    WP_w0(1,var)             = Cond_w(1,var)*n_turns(1,var); % massimo ingombro toroidale WP
end
%% Cechck geometrico sulle dimensioni ottenute nel cavo
r_cable(1:n_layers) = Cond_w(1:n_layers)./Cond_h(1:n_layers);             % Aspect ratio cavi di ogni grades
%% Dati WP    
WP_w = WP_w0(1,1);
WP_w_tot = WP_w+2*GoundIns;
WP_h = sum(Cond_h(1:n_layers));
WP_h_tot = WP_h+2*GoundIns;
A_WP = sum(Cond_h(1:n_layers).*Cond_w(1:n_layers).*n_turns);
A_CICC_tot = sum(S_CICC(1:n_layers).*n_turns);  % Sezione totale cavo non isolato
A_SC_tot = sum(S_Cable(1:n_layers).*n_turns);   % Sezione totale S/C+stabilizer 
A_JT_tot = sum(S_JT(1:n_layers).*n_turns);      % Sezione totale acciaio nei jacket  
Ri_ = R_TF_Innerleg;                            % Quota a raggio esterno Case
Rj_ = Ri_-WP_h-dr_plasma_side-GoundIns*2;       % Quota a raggio interno WP con Ground     

%% Primary radial stress (Pm+Pb)
p_rs = B_TF^2/(2*Mu_0);  % Magnetic pressure, thin WP         
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
r_steel = (Cond_w(n_layers))/(2*JT(n_layers));
S_rm = p_rs*r_steel*scf*dcr_jckt; % Radial memebrane stress Jacket innermost layer  

%% Cechk membrane stress on Jacket 
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
    S_c_VT = 2/(1-beta^2)*p_rs*dcr_vault_tor*1.1; % @ Rk_ correggere con 

    k_bf = 0.5*log(RTFo(dp)/RTFi(dp)); % k bending free
    T_bf = 0.5*(k_bf*n_TF*(n_spire_(1,1)*Iop)^2*Mu_0/(2*pi)); % Hoop tension along TF longitudinal axis

    S_z =  T_bf/(A_JT_tot+A_CASE);  

    S_T_VT = (S_z+S_c_VT); % Vault Tresca stress 
    S_T_JT = (S_z+S_rm_JT);  % Jacket Tresca stress      
end   
V = L*Iop/Tau_discharge*1e-3; % Tensione max singolo TF

Rk_ = 2.88;

%%
% figure('units','normalized','outerposition',[0.25 0.1 0.5 0.9])
figure('units','normalized','outerposition',[0 0 1 1])
% c = uisetcolor([1 1 0],'Select a color')
map = jet(n_layers); 
J_SC = Iop./S_Cable(1:n_layers).*1e-6;
fill([-CASE_w_l/2 -CASE_w/2 CASE_w/2 CASE_w_l/2 -CASE_w_l/2],[Rk_ Ri_ Ri_ Rk_ Rk_],[0.7 0.7 0.7]), hold on
for var = 1:n_layers
    for t = 1:n_turns(var)                                
        % map(var,:) = [min(J_SC)/J_SC(var)    min(J_SC)/J_SC(var)^2    min(J_SC)/J_SC(var)];  
        rectangle('Position',[(-Cond_w(1,var)-Cond_w(1,var)*n_turns(var)/2)+Cond_w(1,var)*t Ri(1,var) Cond_w(1,var) Cond_h(1,var)],...
            'FaceColor',[0.0588    1.0000    1.0000])
        hold on
        rectangle('Position',[(-Cond_w(1,var)-Cond_w(1,var)*n_turns(var)/2+tins(1,var))+Cond_w(1,var)*t Ri(1,var)+tins(1,var) Cond_w(1,var)-2*tins(1,var) Cond_h(1,var)-2*tins(1,var)],...
            'FaceColor',[0.5    0.5    0.5])
        hold on                     
        if shape_cable == 200                  
            rectangle('Position',[(Cond_w(1,var)*(-1-n_turns(var)/2)+Cond_w(1,var)/2-SC_w(1,var)/2)+Cond_w(1,var)*t Ri(1,var)+JT(1,var)+2*tins(1,var)-2*INS_grades SC_w(1,var) SC_w(1,var)],...
            'FaceColor',map(n_layers-var+1,:),'Curvature',[1 1])          
        else
            rectangle('Position',[(-Cond_w(1,var)-Cond_w(1,var)*n_turns(var)/2+JT(1,var)+tins(1,var))+Cond_w(1,var)*t Ri(1,var)+JT(1,var)+tins(1,var) Cond_w(1,var)-2*JT(1,var)-2*tins(1,var) Cond_h(1,var)-2*JT(1,var)-2*tins(1,var)],...
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

TitleName2 = sprintf('n grades=%G, N=%G , Iop=%G kA, Sy(VT)=%G MPa, Sy(JT)=%G MPa, tau=%G s',n_grades,n_spire_(1), round(Iop*1e-3),round(S_amm_VT*1e-6),round(S_amm_JT*1e-6),round(Tau_discharge));
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
text(-dx_*1.03,Rj_+WP_h/2,sprintf('WP(h)=%G m',round(WP_h+GoundIns*2,2)),'fontsize',12,'Interpreter','latex')
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
J_CICC = Iop./(S_Cable(1:n_layers)).*1e-6;
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
            rectangle('Position',[(Cond_w(1,var)*(-1-n_turns(var)/2)+tins(1,var))+Cond_w(1,var)*t Ri(1,var)+tins(1,var) Cond_w(1,var)-2*tins(1,var) Cond_h(1,var)-2*tins(1,var)],...
                'FaceColor',map(colorvar(varg),:))
           hold on
            if shape_cable == 200                   
                rectangle('Position',[(Cond_w(1,var)*(-1-n_turns(var)/2)+Cond_w(1,var)/2-SC_w(1,var)/2)+Cond_w(1,var)*t Ri(1,var)+JT(1,var)+2*tins(1,var)-2*INS_grades SC_w(1,var) SC_w(1,var)],...
                'FaceColor',map(colorvar(varg),:),'Curvature',[1 1])          
            else
                rectangle('Position',[(-Cond_w(1,var)-Cond_w(1,var)*n_turns(var)/2+JT(1,var)+tins(1,var))+Cond_w(1,var)*t Ri(1,var)+JT(1,var)+tins(1,var) Cond_w(1,var)-2*JT(1,var)-2*tins(1,var) Cond_h(1,var)-2*JT(1,var)-2*tins(1,var)],...
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
TitleName1 = sprintf('n grades=%G, N=%G , Iop=%G kA, Sy(VT)=%G MPa, Sy(JT)=%G MPa, tau=%G s ',n_grades,n_spire_(1), round(Iop*1e-3),round(S_amm_VT*1e-6),round(S_amm_JT*1e-6),round(Tau_discharge));
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
text(-dx_*1.03,Rj_+WP_h/2,sprintf('WP(h)=%G m',round(WP_h+GoundIns*2,2)),'fontsize',12,'Interpreter','latex')
text(-dx_*1.03,(Rk_+Rj_)/2,sprintf('Nose(h)=%G m',round(Rj_-Rk_,2)),'fontsize',12,'Interpreter','latex')
text(-dx_*1.03,((Rj_+WP_h)+Ri_)/2,sprintf('dps(h)=%G m',round(dr_plasma_side,2)),'fontsize',12,'Interpreter','latex')
text(0,Ri_*1.1,sprintf('Case(w)=%G m',round(CASE_w,2)),'HorizontalAlignment','center','fontsize',12,'Interpreter','latex')
text(0,Ri_*1.05,sprintf('WP(w)=%G m',round(WP_w,2)),'HorizontalAlignment','center','fontsize',12,'Interpreter','latex')

if printcond == 1, print('-dpng','-r300',TitleName1), end
%% Plot WP - Jeng comparison
% 
% figure('units','normalized','outerposition',[0 0 1 1])
% jeng_ = ceil(Iop./(Cond_w(1).*Cond_h(1)).*1e-6);
% J_CICC = Iop./(S_Cable(1:n_layers)).*1e-6;
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
%         rectangle('Position',[(Cond_w(1,var)*(-1-n_turns(var)/2)+tins(1,var))+Cond_w(1,var)*t Ri(1,var)+tins(1,var) Cond_w(1,var)-2*tins(1,var) Cond_h(1,var)-2*tins(1,var)],...
%             'FaceColor',map(jeng_-mod_min,:))
%        hold on
%         if shape_cable == 200                      
%             rectangle('Position',[(Cond_w(1,var)*(-1-n_turns(var)/2)+Cond_w(1,var)/2-SC_w(1,var)/2)+Cond_w(1,var)*t Ri(1,var)+JT(1,var)+2*tins(1,var)-2*INS_grades SC_w(1,var) SC_w(1,var)],...
%             'FaceColor',map(jeng_-mod_min,:),'Curvature',[1 1])          
%         else
%             rectangle('Position',[(-Cond_w(1,var)-Cond_w(1,var)*n_turns(var)/2+JT(1,var)+tins(1,var))+Cond_w(1,var)*t Ri(1,var)+JT(1,var)+tins(1,var) Cond_w(1,var)-2*JT(1,var)-2*tins(1,var) Cond_h(1,var)-2*JT(1,var)-2*tins(1,var)],...
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
% %% Plot WP - B comparison
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
%         rectangle('Position',[(Cond_w(1,var)*(-1-n_turns(var)/2)+tins(1,var))+Cond_w(1,var)*t Ri(1,var)+tins(1,var) Cond_w(1,var)-2*tins(1,var) Cond_h(1,var)-2*tins(1,var)],...
%             'FaceColor',[0.5    0.5    0.5])
%        hold on
%         if shape_cable == 200                      
%             rectangle('Position',[(Cond_w(1,var)*(-1-n_turns(var)/2)+Cond_w(1,var)/2-SC_w(1,var)/2)+Cond_w(1,var)*t Ri(1,var)+JT(1,var)+2*tins(1,var)-2*INS_grades SC_w(1,var) SC_w(1,var)],...
%             'FaceColor',map(jeng_-mod_min,:),'Curvature',[1 1])          
%         else
%             rectangle('Position',[(-Cond_w(1,var)-Cond_w(1,var)*n_turns(var)/2+JT(1,var)+tins(1,var))+Cond_w(1,var)*t Ri(1,var)+JT(1,var)+tins(1,var) Cond_w(1,var)-2*JT(1,var)-2*tins(1,var) Cond_h(1,var)-2*JT(1,var)-2*tins(1,var)],...
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


