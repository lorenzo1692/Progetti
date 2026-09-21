%% validate_TF_FEM_benchmark_2026.m
% Forward-evaluation check of the analytical Jacket/Case stress formulas
% against the ANSYS APDL FEM benchmark (TF_FEM_benchmark_v0_1, run of
% 21-Sep-2026, patched geometry: fixed H_GRADES indexing and Rj_ closure).
%
% This is a POINT evaluation, not a call into search/scan_wp_designs.m:
% the FEM geometry (100 turns, 2 grades, 11 rows, JT=3.5mm, DTF=135.4mm)
% is fixed and known, so instead of letting size_cicc_cable.m/
% size_case_vault.m iterate to size JT/DTF, this script evaluates the same
% underlying formulas ONCE at that fixed geometry and compares the
% resulting Jacket/Case stress against the confirmed FEM peaks. Keep the
% formulas below in sync with physics/size_cicc_cable.m (cavo_stiffness
% subfunction) and physics/size_case_vault.m if those change.
%
% Confirmed FEM reference (TFBM_jacket.txt / TFBM_case.txt, nodal-averaged,
% global Cartesian axes, LOAD STEP=1 SUBSTEP=1):
%   Jacket (mat 2), node 183613: SINT = 979.6 MPa, SEQV = 898.0 MPa
%   Case   (mat 10), node 217646: SINT = 773.6 MPa, SEQV = 698.4 MPa
% (discretization-error bounds: Jacket up to ~1095 MPa, Case up to ~792 MPa)

clearvars; clc

%% FEM-confirmed geometry and loads (do not edit without a new export)
n_TF = 12;
theta_TF = 2*pi/n_TF;
Mu_0 = 4e-7*pi;

Ri_ = 1.18;                 % Case outer radius (RI_)
dr_plasma_side = 0.02;      % CASE_THICK
GoundIns = 0.005;           % GIT
CASE_w = 2*Ri_*tan(theta_TF/2);

n_turns = [10 10 10 10 10 10 8 8 8 8 8];   % 6 rows grade1 (10 t/row) + 5 rows grade2 (8 t/row)
Cond_h  = [repmat(0.0326,1,6) repmat(0.0228,1,5)];
Cond_w  = repmat(0.046,1,11);
JT      = repmat(0.0035,1,11);
tins    = repmat(0.001,1,11);
n_layers = 11;

E_jckt = 205e9;
E_case = 205e9;
E_ins  = 12e9;

B_TF = 13.489;              % BSUM peak on conductor (TFBM_em.txt / EM_2D009.png)
Iop  = 62300;
n_spire1 = sum(n_turns);    % 100
p_rs = B_TF^2/(2*Mu_0);

Rj_fem = 0.8354;            % corrected Rj_ (includes CASE_THICK, see patches/change_notes.json)
Rk_fem = 0.70;               % RK_
DTF_fem = Rj_fem - Rk_fem;

T_bf = 32675062.45;         % GSBDATA axial force (CHECK_TF.txt / benchmark_case.json)

Re1 = Ri_ - dr_plasma_side - GoundIns;
CASE_w_at_Re1 = 2*Re1*tan(theta_TF/2);
lateral_w = (CASE_w_at_Re1 - Cond_w(1)*n_turns(1) - 2*GoundIns)/2;

%% Reference FEM stress (confirmed, see header)
FEM_Jacket_SINT = 979.6e6; FEM_Jacket_SEQV = 898.0e6;
FEM_Case_SINT   = 773.6e6; FEM_Case_SEQV   = 698.4e6;

fprintf('%-55s %10s %10s\n','Test','S_T_JT','S_T_VT');
run_case(10e9, 0.004, 'A: FEM-matched (E_cbl=10GPa, r_SC=4mm)', 'LTS', ...
    n_layers,Cond_h,Cond_w,JT,tins,n_turns,E_jckt,E_case,E_ins,p_rs,T_bf, ...
    Ri_,CASE_w,theta_TF,lateral_w,Rj_fem,Rk_fem,DTF_fem);
run_case(0.1e9, 0.005, 'B: current tool defaults (E_cbl=0.1GPa, r_SC=5mm)', 'LTS', ...
    n_layers,Cond_h,Cond_w,JT,tins,n_turns,E_jckt,E_case,E_ins,p_rs,T_bf, ...
    Ri_,CASE_w,theta_TF,lateral_w,Rj_fem,Rk_fem,DTF_fem);

fprintf('\nFEM reference: Jacket SINT=%.1f MPa SEQV=%.1f MPa | Case SINT=%.1f MPa SEQV=%.1f MPa\n', ...
    FEM_Jacket_SINT/1e6, FEM_Jacket_SEQV/1e6, FEM_Case_SINT/1e6, FEM_Case_SEQV/1e6);

function run_case(E_cbl, r_SC, label, type_cable, n_layers,Cond_h,Cond_w,JT,tins,n_turns, ...
    E_jckt,E_case,E_ins,p_rs,T_bf,Ri_,CASE_w,theta_TF,lateral_w,Rj_,Rk_,DTF)

SC_h = Cond_h - 2*JT - 2*tins;
SC_w = Cond_w - 2*JT - 2*tins;
R_J = r_SC + JT;
S_Cable = SC_h.*SC_w - (4-pi)*r_SC^2;
S_CICC = (Cond_w-2*tins).*(Cond_h-2*tins) - (4-pi)*R_J.^2;
S_JT = S_CICC - S_Cable;

Ke_rad = zeros(1,n_layers); Ke_tor = zeros(1,n_layers);
for i = 1:n_layers
    [Ke_rad(i), Ke_tor(i)] = cavo_stiffness(E_jckt, E_cbl, E_ins, JT(i), tins(i), Cond_h(i), Cond_w(i), SC_w(i), SC_h(i));
end

WP_h = sum(Cond_h);
A_WP = sum(Cond_h.*Cond_w.*n_turns);
A_SC_tot = sum(S_Cable.*n_turns);
A_JT_tot = sum(S_JT.*n_turns);

%% Primary radial stress (Pm+Pb) - evaluated at every layer, keeping the
% worst case (matches search/scan_wp_designs.m). Layers 6/7 (the grade1/
% grade2 row boundary) get the provisional transition SCF - see
% TF_FEM_benchmark_2026_findings.md for the calibration and its caveats.
SCF_transition_provisional = 3.15;
is_transition = false(1,n_layers);
is_transition([6 7]) = true;

xxx = [1.087,1.136,1.190,1.250,1.316,1.389,1.471,1.563,1.667,1.786,1.923,2.083,2.273,2.500,2.778];
if strcmp(type_cable,'LTS')
    yyy = [1.01,1.03,1.06,1.10,1.16,1.22,1.29,1.36,1.43,1.50,1.57,1.64,1.71,1.78,1.85];
else
    yyy = [1.01,1.02,1.02,1.04,1.06,1.07,1.11,1.13,1.16,1.18,1.22,1.27,1.29,1.29,1.31];
end
S_rm_per_layer = zeros(1,n_layers);
for var = 1:n_layers
    param = Cond_w(var)/SC_w(var);
    if param > 1 && param < 2.8
        pf = polyfit(xxx,yyy,5);
        scf = polyval(pf,param);
    else
        scf = 1.5;
    end
    K_jckt = 2*JT(var)/Cond_h(var)*E_jckt;
    dcr_jckt = K_jckt/Ke_rad(var);
    r_steel = (Cond_w(var)-2*tins(var))/(2*JT(var));
    S_rm_var = p_rs*r_steel*scf*dcr_jckt;
    if is_transition(var)
        S_rm_var = S_rm_var*SCF_transition_provisional;
    end
    S_rm_per_layer(var) = S_rm_var;
end
S_rm = max(S_rm_per_layer);

%% Case/vault, single evaluation at the FEM-confirmed DTF (see size_case_vault.m)
CASE_w_l = 2*Rk_*tan(theta_TF/2);
A_tot = (CASE_w + CASE_w_l)*(Ri_ - Rk_)/2;
A_CASE = A_tot - A_WP;

Ke_WP_rad = sum(1./(Ke_rad.*n_turns))^-1;
K_ps_rad = E_case/0.02*(2*Ri_*tan(theta_TF/2)); % dr_plasma_side hardcoded (CASE_THICK)
K_lat_rad = E_case*(lateral_w/2)/WP_h;
K_vault_rad = E_case*CASE_w_l/DTF;
Ke_case_rad = (1/(Ke_WP_rad+2*K_lat_rad) + 1/K_vault_rad + 1/K_ps_rad)^-1;
dcr_WP_rad = Ke_WP_rad/Ke_case_rad;
S_rm_JT = S_rm*dcr_WP_rad;

h_unit = 1;
k_steel_tor = E_case*(h_unit*(Ri_-Rk_)*2*pi*(Ri_+Rk_)/2);
k_SC_tor = E_cbl*(h_unit*(Ri_-Rj_)*2*pi*(Ri_+Rj_)/2)*(1-(A_WP-A_SC_tot)/A_WP);
k_JT_tor = E_jckt*(h_unit*(Ri_-Rj_)*2*pi*(Ri_+Rj_)/2)*(1-(A_WP-A_JT_tot)/A_WP);
k_vault_tor = E_case*(h_unit*(Rj_-Rk_)*2*pi*(Rj_+Rk_)/2);
dcr_vault_tor = k_steel_tor/(k_SC_tor+k_JT_tor+k_vault_tor);

beta = Rk_/Ri_;
S_c_VT = 2/(1-beta^2)*p_rs*dcr_vault_tor;
S_z = T_bf/(A_JT_tot+A_CASE);

S_T_VT = S_z + S_c_VT;
S_T_JT = S_z + S_rm_JT;

fprintf('%-55s %8.1f M %8.1f M\n', label, S_T_JT/1e6, S_T_VT/1e6);
end

function [Ke_rad, Ke_tor] = cavo_stiffness(E_jckt, E_cbl, E_ins, JT, tins, Cond_h, Cond_w, SC_w, SC_h)
Ke_rad = 2*E_jckt*JT/Cond_h + 2*tins*E_ins/Cond_h + ...
    (1/(E_cbl*SC_w/SC_h) + 2/(E_jckt*Cond_w/JT) + 2/(E_ins*Cond_w/tins))^-1;
Ke_tor = 2*E_jckt*JT/Cond_w + 2*tins*E_ins/Cond_w + ...
    (1/(E_cbl*SC_h/SC_w) + 2/(E_jckt*JT/Cond_w) + 2/(E_ins*tins/Cond_w))^-1;
end
