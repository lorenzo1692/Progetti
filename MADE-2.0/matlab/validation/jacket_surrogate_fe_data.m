function d = jacket_surrogate_fe_data(row, p, name, opts)
%JACKET_SURROGATE_FE_DATA FE jacket stresses per layer + the scan's analytic quantities.
%   d = JACKET_SURROGATE_FE_DATA(row, p, name) runs WP_MECH_SURROGATE on a
%   design point (primary and total load cases) and returns, per layer, the
%   jacket Pm and Pm+Pb together with the quantities the scan uses for the
%   jacket (sigma_nom = p_rs*r_steel*dcr_jckt, S_z, dcr_WP_rad, layer field,
%   turns): the data behind JACKET_STRESS_SURROGATE (see
%   validation/tools/fit_jacket_surrogate.py and
%   validation/results/jacket_surrogate_calibration.csv, one row per layer).
if nargin < 4, opts = struct('verbose', false); end
Mu_0 = 4e-7*pi;
g = compute_operating_params(p);
nl = row.n_layers; nt = row.n_turns(1:nl);
tg = wp_turn_geometry(row, p);
tins = tg.tins;
Cw = row.Cond_w(1:nl); Ch = row.Cond_h(1:nl); JT = row.JT(1:nl);
% cable type / E per layer
E_cbl = p.E_cbl_LTS*ones(1, nl);
if isfield(row, 'type_cable')
    tc = row.type_cable; for k = 1:nl, if strcmp(tc{k}, 'HTS'), E_cbl(k) = p.E_cbl_HTS; end, end
end
SCw = tg.cab_w; SCh = tg.cab_h;
Ke = 2*p.E_jckt*JT./Ch + 2*tins*p.E_ins./Ch + 1./(1./(E_cbl.*SCw./SCh) + 2./(p.E_jckt*Cw./JT) + 2./(p.E_ins*Cw./tins));
K_jckt = 2*p.E_jckt*JT./Ch;
dcr_jckt = K_jckt./Ke;
r_steel = (Cw - 2*tins)./(2*JT);
param = Cw./SCw;
% layer peak field (discrete) and pressure as the scan
Bl = wp_peak_field_fast(row, p);
if isfield(row, 'B_grade') && all(isfinite(row.B_grade(1))), Bp = row.B_grade(1); else, Bp = Bl(1); end
p_rs = Bp^2/(2*Mu_0);
sigma_nom = p_rs*r_steel.*dcr_jckt;
% transitions: cable changes between adjacent layers
chg = abs(diff(Ch)) > 1e-9 | abs(diff(JT)) > 1e-9;
if isfield(row, 'S_Cable'), sc = row.S_Cable(1:nl); chg = chg | abs(diff(sc)) > 1e-12; end
T = false(1, nl); idx = find(chg); T(idx) = true; T(idx+1) = true;
% accumulated radial load fraction (Lorentz force ~ n_k*B_k)
f = nt.*Bl; A = (cumsum(f) - 0.5*f)/sum(f);
% vault / axial quantities of size_case_vault at the final nose
theta = 2*pi/p.n_TF; Ri_ = row.Ri_; Rk_ = row.Rk_;
WP_h = sum(Ch); Rj_ = Ri_ - WP_h - (nl-1)*p.INS_grades - p.dr_plasma_side - 2*p.GoundIns;
DTF = Rj_ - Rk_;
CASE_w = 2*Ri_*tan(theta/2);
CASE_w_l = 2*Rk_*tan(theta/2);
A_WP = sum(Ch.*Cw.*nt); A_JT_tot = sum(tg.A_jacket.*nt);
A_CASE = (CASE_w + CASE_w_l)*(Ri_ - Rk_)/2 - A_WP;
Re1 = Ri_ - p.dr_plasma_side - p.GoundIns;
lateral_w = (2*Re1*tan(theta/2) - 2*p.GoundIns - Cw(1)*nt(1))/2;   % as the scan's WP_w0
if isfield(row, 'lateral_w'), lateral_w = row.lateral_w; end
Ke_WP_rad = sum(1./(Ke.*nt))^-1;
K_ps_rad = p.E_case/p.dr_plasma_side*(2*Ri_*tan(theta/2));
K_lat_rad = p.E_case*(lateral_w/2)/WP_h;
K_vault_rad = p.E_case*CASE_w_l/DTF;
Ke_case_rad = (1/(Ke_WP_rad + 2*K_lat_rad) + 1/K_vault_rad + 1/K_ps_rad)^-1;
dcr_WP_rad = Ke_WP_rad/Ke_case_rad;
T_bf = 0.5*(g.k_bf*p.n_TF*(sum(nt)*row.Iop)^2*Mu_0/(2*pi));
S_z = T_bf/(A_JT_tot + A_CASE);
% FE
opts.classify = 1;
t0 = tic; out = wp_mech_surrogate(row, p, opts);
PmPb_P = out.primary.layer.PmPb; PmPb_T = out.layer.PmPb; Pm_P = out.primary.layer.Pm; Pm_T = out.layer.Pm;
SCF_FE = (PmPb_P - S_z)./(dcr_WP_rad*sigma_nom);
d = struct('name', name, 'n_layers', nl, 'n_turns', nt, 'W', max(nt.*Cw), 'Iop', row.Iop, ...
    'param', param, 'A', A, 'T', T, 'Bl', Bl, 'p_rs', p_rs, 'sigma_nom', sigma_nom, ...
    'dcr_jckt', dcr_jckt, 'r_steel', r_steel, 'JT_Ch', JT./Ch, 'Ch_Cw', Ch./Cw, 'S_z', S_z, 'dcr_WP_rad', dcr_WP_rad, ...
    'PmPb_P', PmPb_P, 'PmPb_T', PmPb_T, 'Pm_P', Pm_P, 'Pm_T', Pm_T, 'SCF_FE', SCF_FE, 'valid', out.valid, 'checks', out.checks.summary, ...
    'fe_time', toc(t0));
fprintf('%s: %d layers, valid %d (%s), FE %.0f s, SCF_FE %s\n', name, nl, out.valid, out.checks.summary, ...
    d.fe_time, mat2str(round(SCF_FE*100)/100));
end
