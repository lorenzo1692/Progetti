function out = forward_eval_wp_stress(geo, mat, B, T_bf, is_transition, SCF_transition)
%FORWARD_EVAL_WP_STRESS Analytical Jacket/Case stress at a FIXED WP geometry.
%
%   out = FORWARD_EVAL_WP_STRESS(geo, mat, B, T_bf, is_transition, SCF_transition)
%   evaluates, once and with no sizing iteration, the same formulas that
%   physics/size_cicc_cable.m (cable stiffness), search/scan_wp_designs.m
%   (per-layer radial membrane stress) and physics/size_case_vault.m
%   (case/vault stress balance) use, at a geometry that is already known
%   (e.g. the one a FEM run was built on). Used to compare the analytical
%   model against FEM layer by layer. Keep in sync with those files.
%
%   geo - n_turns, Cond_h, Cond_w, JT, tins (per layer), Ri_, Rj_, Rk_,
%         lateral_w, n_TF, dr_plasma_side, r_SC_min, r_SC_max
%   mat - E_jckt, E_case, E_ins, E_cbl                         [Pa]
%   B   - peak field used for the magnetic pressure p = B^2/(2 mu0) [T]
%   T_bf - axial (vertical) force per leg                      [N]
%   is_transition  - logical per layer, layers that get the transition SCF
%   SCF_transition - multiplier applied to those layers (1 = none)
%
%   out - S_T_JT (per layer), S_T_JT_max, S_T_VT, S_z, S_rm_JT (per layer),
%         dcr_WP_rad, DTF                                        [Pa / m]

Mu_0 = 4e-7*pi;
theta_TF = 2*pi/geo.n_TF;
n_layers = numel(geo.Cond_h);
nt = geo.n_turns; Ch = geo.Cond_h; Cw = geo.Cond_w; JT = geo.JT; tins = geo.tins;

r_SC = min(max(JT, geo.r_SC_min), geo.r_SC_max);
SC_h = Ch - 2*JT - 2*tins;
SC_w = Cw - 2*JT - 2*tins;
R_J = r_SC + JT;
S_Cable = SC_h.*SC_w - (4-pi)*r_SC.^2;
S_CICC = (Cw-2*tins).*(Ch-2*tins) - (4-pi)*R_J.^2;
S_JT = S_CICC - S_Cable;

Ke_rad = 2*mat.E_jckt*JT./Ch + 2*tins*mat.E_ins./Ch + ...
    (1./(mat.E_cbl*SC_w./SC_h) + 2./(mat.E_jckt*Cw./JT) + 2./(mat.E_ins*Cw./tins)).^-1;

WP_h = sum(Ch);
A_WP = sum(Ch.*Cw.*nt);
A_SC_tot = sum(S_Cable.*nt);
A_JT_tot = sum(S_JT.*nt);
p_rs = B^2/(2*Mu_0);

% Per-layer radial membrane stress (scan_wp_designs.m)
xxx = [1.087,1.136,1.190,1.250,1.316,1.389,1.471,1.563,1.667,1.786,1.923,2.083,2.273,2.500,2.778];
yyy = [1.01,1.03,1.06,1.10,1.16,1.22,1.29,1.36,1.43,1.50,1.57,1.64,1.71,1.78,1.85];
pf = polyfit(xxx, yyy, 5);
S_rm = zeros(1, n_layers);
for k = 1:n_layers
    param = Cw(k)/SC_w(k);
    if param > 1 && param < 2.8
        scf = polyval(pf, param);
    else
        scf = 1.5;
    end
    dcr_jckt = (2*JT(k)/Ch(k)*mat.E_jckt)/Ke_rad(k);
    r_steel = (Cw(k)-2*tins(k))/(2*JT(k));
    S_rm(k) = p_rs*r_steel*scf*dcr_jckt;
    if is_transition(k)
        S_rm(k) = S_rm(k)*SCF_transition;
    end
end

% Case/vault (size_case_vault.m, single evaluation at the given DTF)
CASE_w = 2*geo.Ri_*tan(theta_TF/2);
CASE_w_l = 2*geo.Rk_*tan(theta_TF/2);
DTF = geo.Rj_ - geo.Rk_;
A_CASE = (CASE_w + CASE_w_l)*(geo.Ri_ - geo.Rk_)/2 - A_WP;
Ke_WP_rad = sum(1./(Ke_rad.*nt))^-1;
K_ps_rad = mat.E_case/geo.dr_plasma_side*CASE_w;
K_lat_rad = mat.E_case*(geo.lateral_w/2)/WP_h;
K_vault_rad = mat.E_case*CASE_w_l/DTF;
Ke_case_rad = (1/(Ke_WP_rad+2*K_lat_rad) + 1/K_vault_rad + 1/K_ps_rad)^-1;
dcr_WP_rad = Ke_WP_rad/Ke_case_rad;

Ri_ = geo.Ri_; Rj_ = geo.Rj_; Rk_ = geo.Rk_;
k_steel_tor = mat.E_case*((Ri_-Rk_)*2*pi*(Ri_+Rk_)/2);
k_SC_tor = mat.E_cbl*((Ri_-Rj_)*2*pi*(Ri_+Rj_)/2)*(1-(A_WP-A_SC_tot)/A_WP);
k_JT_tor = mat.E_jckt*((Ri_-Rj_)*2*pi*(Ri_+Rj_)/2)*(1-(A_WP-A_JT_tot)/A_WP);
k_vault_tor = mat.E_case*((Rj_-Rk_)*2*pi*(Rj_+Rk_)/2);
dcr_vault_tor = k_steel_tor/(k_SC_tor+k_JT_tor+k_vault_tor);
beta = Rk_/Ri_;
S_c_VT = 2/(1-beta^2)*p_rs*dcr_vault_tor;
S_z = T_bf/(A_JT_tot + A_CASE);

out.S_rm_JT = S_rm*dcr_WP_rad;
out.S_T_JT = S_z + out.S_rm_JT;
out.S_T_JT_max = max(out.S_T_JT);
out.S_T_VT = S_z + S_c_VT;
out.S_z = S_z;
out.dcr_WP_rad = dcr_WP_rad;
out.DTF = DTF;
end
