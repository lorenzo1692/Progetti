%% verify_made_vs_gps2d.m
% The quantities the SCAN of MADE computes vs the two ANSYS 2D
% generalized-plane-strain runs of 2026 (TF_FEM_benchmark_2026, design 7
% EM_2D007), at the geometry the FEM was built on (no sizing iteration):
%
%   1. peak field on the conductor: smeared Ampere (field_model 0) and the
%      discrete 2D model that the calibrated field (field_model 1) is fitted
%      on (within ~0.1 T, see physics/wp_field_calibration.m) vs FEM BSUM;
%   2. jacket, per layer: the fast surrogate (scf_model 1,
%      physics/jacket_stress_surrogate.m) and the analytic formula with the
%      SCF table and SCF_transition_provisional (scf_model 0) vs the ANSYS
%      linearized Pm and Pm+Pb (same sections as WP_MECH_SURROGATE);
%   3. case: the vault formula of SIZE_CASE_VAULT (S_z + S_c_VT) vs the
%      ANSYS case SCLs.
%
% LOAD CASES. The two ANSYS runs apply cool-down + Lorentz + axial force in
% one step: their linearized stresses are TOTAL (P+Q). The criteria of the
% scan (Pm <= Sm, Pm+Pb <= 1.5 Sm) and the jacket surrogate are PRIMARY
% (Lorentz + axial, no cool-down). The chain is therefore shown in three
% links: ANSYS total vs 2D FE (WP_MECH_SURROGATE) total, 2D FE primary vs
% the jacket surrogate; the 2D FE values are read from
% validation/results/jacket_surrogate_calibration.csv (cases 'bench', 'd7').
% A direct ANSYS check of the primary stresses needs a run without the
% cool-down (see docs/RIS_E_RESIDUO_FEM.md).
%
% The axial force T_bf is taken from each FEM run (like with like). The
% 2D FE surrogate itself is compared with the same ANSYS data by
% validate_mech_surrogate_2026.m. Note: the jacket surrogate coefficients
% were fitted on WP_MECH_SURROGATE results that include these two
% geometries, so for the jacket this is a check of the chain
% scan -> 2D FE -> ANSYS on the calibration points, not an independent test.

clearvars; clc
this_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(this_dir, '..')); made_paths('TF');
Mu_0 = 4e-7*pi;
pp = struct('n_TF', 12, 'dr_plasma_side', 0.02, 'GoundIns', 0.005, 'INS_grades', 0.0005, ...
    'turn_insulation_nominal', 0.001, 'Increm', 1, 'shape_cable', 201);   % both FEM runs: Rect
mat = struct('E_jckt', 205e9, 'E_case', 205e9, 'E_ins', 12e9, 'E_cbl', 10e9);
SCF_tr = 3.15;
scf_x = [1.087,1.136,1.190,1.250,1.316,1.389,1.471,1.563,1.667,1.786,1.923,2.083,2.273,2.500,2.778];
scf_y = [1.01,1.03,1.06,1.10,1.16,1.22,1.29,1.36,1.43,1.50,1.57,1.64,1.71,1.78,1.85];

cases = struct([]);
% --- A. benchmark (ANSYS 21-Sep-2026) --------------------------------------
c = struct();
c.name = 'TF_FEM_benchmark_2026'; c.tag = 'bench';
c.row = struct('Iop', 62300, 'n_layers', 11, 'n_turns', [10 10 10 10 10 10 8 8 8 8 8], ...
    'Cond_w', repmat(0.046,1,11), 'Cond_h', [repmat(0.0326,1,6) repmat(0.0228,1,5)], ...
    'JT', repmat(0.0035,1,11), 'Ri_', 1.18, 'Rk_', 0.70);
c.Rj_ = 0.8354; c.T_bf = 32675062.45; c.grade_first = [1 7]; c.r_SC = [0.004 0.004];
c.B_fem = 13.489;
c.Pm   = [487 472 485 503 526 581 505 516 530 543 549];
c.PmPb = [632 640 648 693 739 863 634 677 717 754 778];
c.scl  = [697 736; 689 700; 670 779; 418 420; 411 413; 389 400; 492 508];
c.case_sint = 773.6;
cases = [cases, c];
% --- B. design 7 (EM_2D007, 24-Sep-2026) ----------------------------------
c = struct();
c.name = 'design 7 (EM_2D007)'; c.tag = 'd7';
c.row = struct('Iop', 64628, 'n_layers', 13, 'n_turns', 8*ones(1,13), ...
    'Cond_w', repmat(0.04300426499905924,1,13), ...
    'Cond_h', [repmat(0.034167368089179834,1,4) repmat(0.026104457386404455,1,3) repmat(0.024226622311504957,1,6)], ...
    'JT', repmat(0.003,1,13), 'Ri_', 1.259425287356322, 'Rk_', 0.7300827089713595);
c.Rj_ = c.row.Ri_ - sum(c.row.Cond_h) - 12*pp.INS_grades - pp.dr_plasma_side - 2*pp.GoundIns;
c.T_bf = 36.266e6; c.grade_first = [1 5 8]; c.r_SC = [0.002 0.006];
c.B_fem = 14.482;
c.Pm   = [474 459 481 508 503 524 540 543 552 560 568 573 581];
c.PmPb = [641 718 791 840 772 815 841 827 844 856 868 889 963];
c.scl  = [639 745; 639 657; 647 789; 363 370; 347 354; 335 339; 505 532];
c.case_sint = 800.5;
cases = [cases, c];
% 2D FE results (WP_MECH_SURROGATE) of the calibration data set
fid = fopen(fullfile(this_dir, 'results', 'jacket_surrogate_calibration.csv'));
hdr = strsplit(fgetl(fid), ',');
C = textscan(fid, ['%s' repmat('%f', 1, numel(hdr)-1)], 'Delimiter', ',');
fclose(fid);
col = @(name) C{strcmp(hdr, name)};
scl_names = {'nose centreline', 'nose mid', 'vault diagonal', 'side wall 25%', 'side wall 50%', ...
    'side wall 75%', 'plasma-side plate'};

for i = 1:numel(cases)
    c = cases(i); r = c.row; nl = r.n_layers;
    fprintf('\n===== %s =====\n', c.name);
    % 1. field
    theta = 2*pi/pp.n_TF;
    Re1 = r.Ri_ - pp.dr_plasma_side - pp.GoundIns;
    B_amp = Mu_0*pp.n_TF*sum(r.n_turns)*r.Iop/(2*pi*Re1);   % smeared, at the plasma-side WP edge
    r.B_TF = B_amp;
    fld = compute_discrete_field_profile(r, pp);
    B_lay = arrayfun(@(k) max(fld.B_discrete(fld.layer == k)), 1:nl);
    fprintf('peak field on conductor: FEM %.3f T | discrete (= calibrated) %.3f T (%+.1f %%) | smeared Ampere %.3f T (%+.1f %%)\n', ...
        c.B_fem, max(B_lay), 100*(max(B_lay)/c.B_fem - 1), B_amp, 100*(B_amp/c.B_fem - 1));
    % 2. jacket - scan quantities at the FEM geometry
    tins = pp.turn_insulation_nominal*ones(1, nl);
    lateral_w = (2*Re1*tan(theta/2) - r.Cond_w(1)*r.n_turns(1) - 2*pp.GoundIns)/2;
    geo = struct('n_turns', r.n_turns, 'Cond_h', r.Cond_h, 'Cond_w', r.Cond_w, 'JT', r.JT, 'tins', tins, ...
        'Ri_', r.Ri_, 'Rj_', c.Rj_, 'Rk_', r.Rk_, 'lateral_w', lateral_w, 'n_TF', pp.n_TF, ...
        'dr_plasma_side', pp.dr_plasma_side, 'r_SC_min', c.r_SC(1), 'r_SC_max', c.r_SC(2));
    is_tr = false(1, nl);
    for k = 2:numel(c.grade_first), is_tr([c.grade_first(k)-1 c.grade_first(k)]) = true; end
    B1 = max(B_lay);
    o0 = forward_eval_wp_stress(geo, mat, B1, c.T_bf, is_tr, SCF_tr);      % scf_model 0
    o1 = forward_eval_wp_stress(geo, mat, B1, c.T_bf, false(1,nl), 1);     % for sigma_nom, S_z
    p_rs = B1^2/(2*Mu_0);
    sigma_nom = o1.S_rm_JT/o1.dcr_WP_rad;          % p_rs r_steel dcr_jckt SCF(table)
    prm = r.Cond_w./(r.Cond_w - 2*r.JT - 2*tins);           % SCF table of SCAN_WP_DESIGNS (LTS)
    SCFtab = 1.5*ones(1, nl); inr = prm > 1 & prm < 2.8;
    SCFtab(inr) = polyval(polyfit(scf_x, scf_y, 5), prm(inr));
    sigma_nom = sigma_nom./SCFtab;                 % without the SCF: the surrogate input
    [Pm_s, PmPb_s] = jacket_stress_surrogate(sigma_nom, p_rs, r.n_turns, B_lay, r.Iop, ...
        r.Cond_w(1)*r.n_turns(1), o1.S_z);
    fprintf('axial stress S_z %.0f MPa (T_bf from the FEM %.2f MN)\n', o1.S_z/1e6, c.T_bf/1e6);
    sel = strcmp(C{1}, c.tag); kk = col('k'); kk = kk(sel);
    FE = zeros(nl, 4); nm = {'Pm_T', 'PmPb_T', 'Pm_P', 'PmPb_P'};
    for q = 1:4, v = col(nm{q}); FE(kk, q) = v(sel)/1e6; end
    fprintf('Pm+Pb per layer [MPa]:    TOTAL (cool-down + Lorentz + axial)  |  PRIMARY (Lorentz + axial)\n');
    fprintf('%-4s %8s %8s %7s   | %8s %9s %7s | %8s\n', 'lay', 'ANSYS', '2D FE', 'FE/AN', '2D FE', 'surrogate', 'sur/FE', 'scf_mod0');
    for k = 1:nl
        fprintf('L%-3d %8.0f %8.0f %+6.0f%%   | %8.0f %9.0f %+6.0f%% | %8.0f\n', k, c.PmPb(k), FE(k,2), ...
            100*(FE(k,2)/c.PmPb(k) - 1), FE(k,4), PmPb_s(k)/1e6, 100*(PmPb_s(k)/1e6/FE(k,4) - 1), o0.S_T_JT(k)/1e6);
    end
    fprintf('design max Pm   : ANSYS tot %.0f | FE tot %.0f (%+.1f %%) | FE prim %.0f | surrogate %.0f (%+.1f %%)\n', ...
        max(c.Pm), max(FE(:,1)), 100*(max(FE(:,1))/max(c.Pm) - 1), max(FE(:,3)), max(Pm_s)/1e6, ...
        100*(max(Pm_s)/1e6/max(FE(:,3)) - 1));
    fprintf('design max Pm+Pb: ANSYS tot %.0f | FE tot %.0f (%+.1f %%) | FE prim %.0f | surrogate %.0f (%+.1f %%) | scf_model 0 %.0f (%+.1f %% vs FE prim)\n', ...
        max(c.PmPb), max(FE(:,2)), 100*(max(FE(:,2))/max(c.PmPb) - 1), max(FE(:,4)), max(PmPb_s)/1e6, ...
        100*(max(PmPb_s)/1e6/max(FE(:,4)) - 1), o0.S_T_JT_max/1e6, 100*(o0.S_T_JT_max/1e6/max(FE(:,4)) - 1));
    fprintf('cool-down effect in the 2D FE (total/primary - 1) on the design max Pm+Pb: %+.1f %%\n', ...
        100*(max(FE(:,2))/max(FE(:,4)) - 1));
    % 3. case
    fprintf('case: vault formula S_z + S_c_VT = %.0f MPa | ANSYS nose/vault SCL Pm max %.0f, Pm+Pb max %.0f, peak SINT %.0f MPa\n', ...
        o1.S_T_VT/1e6, max(c.scl(1:3,1)), max(c.scl(1:3,2)), c.case_sint);
    fprintf('      (vault formula vs SCL Pm max %+.1f %%, vs SCL Pm+Pb max %+.1f %%)\n', ...
        100*(o1.S_T_VT/1e6/max(c.scl(1:3,1)) - 1), 100*(o1.S_T_VT/1e6/max(c.scl(1:3,2)) - 1));
    res(i) = struct('name', c.name, 'B_lay', B_lay, 'B_amp', B_amp, 'Pm', Pm_s/1e6, 'PmPb', PmPb_s/1e6, ...
        'FE', FE, 'S_T_JT0', o0.S_T_JT/1e6, 'S_T_VT', o1.S_T_VT/1e6); %#ok<SAGROW>
end

