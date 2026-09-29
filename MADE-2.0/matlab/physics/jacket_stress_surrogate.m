function [Pm, PmPb] = jacket_stress_surrogate(sigma_nom, p_rs, n_turns, B_layer, Iop, W1, S_z)
%JACKET_STRESS_SURROGATE Fast jacket primary Pm and Pm+Pb per layer (scf_model = 1).
%
%   [Pm, PmPb] = JACKET_STRESS_SURROGATE(sigma_nom, p_rs, n_turns, B_layer,
%   Iop, W1, S_z) estimates, for every layer k of a winding pack, the
%   primary (Lorentz + axial, no cool-down) linearized Tresca stresses of
%   the jacket that the 2D FE model WP_MECH_SURROGATE computes (maximum
%   over the turns of the layer), from quantities the scan already has:
%
%     Pm_k   = S_z + a_m * sigma_nom_k + b_m * sigma_acc_k
%     PmPb_k = S_z + a_b * sigma_nom_k + b_b * sigma_acc_k
%
%     sigma_nom_k = p_rs * r_steel_k * dcr_jckt_k   (the analytic jacket
%                   formula of SCAN_WP_DESIGNS without its SCF)
%     sigma_acc_k = q_k / p_rs * sigma_nom_k,  q_k = (sum_{j<k} n_j Iop B_j
%                   + n_k Iop B_k / 2) / W1: the radial line load of the
%                   layers above, per unit toroidal WP width, at the middle
%                   of layer k (the real pressure accumulated towards the
%                   nose, which the uniform p_rs of the analytic formula
%                   does not see: in the FE the jacket stress grows by
%                   ~70% from the plasma-side layer to the deep layers)
%     S_z         = axial stress T_bf / (A_jacket + A_case) (SIZE_CASE_VAULT)
%
%   The case / WP stiffness share dcr_WP_rad of the analytic chain is NOT
%   applied: with it the FE-implied factor varied by 2x between designs
%   with similar FE stresses.
%
%   Inputs (1 x n_layers unless noted): sigma_nom [Pa]; p_rs [Pa] scalar;
%   n_turns; B_layer peak field per layer [T] (the field the layer is
%   sized for); Iop [A]; W1 toroidal width of the first layer [m]; S_z [Pa].
%
%   CALIBRATION (docs/RIS_E_RESIDUO_FEM.md, validation/tools/
%   fit_jacket_surrogate.py, data validation/results/
%   jacket_surrogate_calibration.csv): least squares on 242 layers of 13
%   rectangular-cable designs computed with WP_MECH_SURROGATE (validated
%   against ANSYS within +-5-7% on the linearized stresses): design 7,
%   TF_FEM_benchmark_2026, scan design 10 (29-Sep-2026) and 10 designs of the calibrated-field scan of
%   the template machine (W 304-344 mm, 14-29 layers, 21-66 kA).
%   Leave-one-design-out error on the design maximum: Pm+Pb -25..+13 %
%   (layer rms 12 %), Pm -16..+12 % (layer rms 7 %); the analytic formula
%   with the SCF table and SCF_transition_provisional: Pm+Pb -62..+6 %.
%   Largest under-prediction on narrow WPs whose deep layers are narrowed
%   (local bending peaks at the width steps). A grade-transition term was
%   not significant (coefficient ~0.03) and is not used; a width-step term
%   (turns lost to the next layer) did not reduce the error either.
%   Design 10 (not in the fit when first checked): Pm+Pb 964 vs FE 1075 MPa
%   (-10 %), peak at the edge turn of the last 14-turn layer over the
%   10-turn one. Not calibrated for
%   RIS cables. A PRELIMINARY screening model: the chosen design must be
%   verified with WP_MECH_SURROGATE.

a_m = 0.711; b_m = 0.244;          % Pm
a_b = 1.285; b_b = 0.601;          % Pm + Pb
F = n_turns(:)'.*Iop.*B_layer(:)';
q = (cumsum(F) - 0.5*F)/W1;
sigma_acc = q/p_rs.*sigma_nom(:)';
Pm   = S_z + a_m*sigma_nom(:)' + b_m*sigma_acc;
PmPb = S_z + a_b*sigma_nom(:)' + b_b*sigma_acc;
end
