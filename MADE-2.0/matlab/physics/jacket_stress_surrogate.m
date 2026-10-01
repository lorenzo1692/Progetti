function [Pm, PmPb] = jacket_stress_surrogate(sigma_nom, p_rs, n_turns, B_layer, Iop, W1, S_z, JT)
%JACKET_STRESS_SURROGATE Fast jacket primary Pm and Pm+Pb per layer (scf_model = 1).
%
%   [Pm, PmPb] = JACKET_STRESS_SURROGATE(sigma_nom, p_rs, n_turns, B_layer,
%   Iop, W1, S_z, JT) estimates, for every layer k of a winding pack, the
%   primary (Lorentz + axial, no cool-down) linearized Tresca stresses of
%   the jacket that the 2D FE model WP_MECH_SURROGATE computes (maximum
%   over the turns of the layer), from quantities the scan already has:
%
%     Pm_k   = S_z + g_k*(a_m*sigma_nom_k + b_m*sigma_acc_k + c_m*sigma_acc_k*s_k)
%     PmPb_k = S_z + g_k*(a_b*sigma_nom_k + b_b*sigma_acc_k + c_b*sigma_acc_k*s_k)
%
%     sigma_nom_k = p_rs * r_steel_k * dcr_jckt_k   (the analytic jacket
%                   formula of SCAN_WP_DESIGNS without its SCF; ~ 1/JT)
%     sigma_acc_k = q_k / p_rs * sigma_nom_k,  q_k = (sum_{j<k} n_j Iop B_j
%                   + n_k Iop B_k / 2) / W1: the radial line load of the
%                   layers above, per unit toroidal WP width, at the middle
%                   of layer k (the real pressure accumulated towards the
%                   nose; in the FE the jacket stress grows by ~70% from the
%                   plasma-side layer to the deep layers)
%     s_k         = (n_k - n_{k+1})/n_k, the width step under layer k: the
%                   edge turns of a layer that overhangs a narrower one
%                   carry a local bending peak (0 for the last layer)
%     g_k         = (JT_k/3 mm)^0.2: sigma_nom alone makes the stress fall
%                   ~1/JT, about 1.5x faster than the FE when the jacket of
%                   a given layout is thickened (+50% steel -> -14% in the FE)
%     S_z         = axial stress T_bf / (A_jacket + A_case) (SIZE_CASE_VAULT)
%
%   Without JT (7 inputs) g_k = 1. The case / WP stiffness share dcr_WP_rad
%   of the analytic chain is not applied: with it the FE-implied factor
%   varied by 2x between designs with similar FE stresses.
%
%   Inputs (1 x n_layers unless noted): sigma_nom [Pa]; p_rs [Pa] scalar;
%   n_turns; B_layer peak field per layer [T] (the field the layer is
%   sized for); Iop [A]; W1 toroidal width of the first layer [m]; S_z
%   [Pa]; JT jacket thickness [m].
%
%   CALIBRATION (docs/DIMENSIONAMENTO_JACKET.md, validation/tools/
%   fit_jacket_surrogate.py, data validation/results/
%   jacket_surrogate_calibration.csv): least squares on 392 layers of 23
%   rectangular-cable FE runs (WP_MECH_SURROGATE, validated against ANSYS
%   within +-5-7% on the linearized stresses) of 13 base designs: design
%   7, TF_FEM_benchmark_2026, scan design 10, 10 designs of the
%   calibrated-field scan of the template machine (W 304-544 mm, 11-29
%   layers, 21-66 kA), plus jacket-thickness variants (+0.5, +1.0 mm) of
%   six of them (validation/run_jacket_jt_calibration.m).
%   Leave-one-base-design-out error (all variants of a design left out
%   together), FE/surrogate on the design maximum: Pm+Pb median 1.03,
%   90th percentile 1.17, max 1.23 (layer rms 10%); Pm median 1.03, 90th
%   percentile 1.13, max 1.15 (layer rms 7%). Response to the jacket
%   thickness of a given layout within ~3% of the FE. The scan applies
%   jacket_margin on top (SIZE_JACKET_SURROGATE). The residual scatter
%   comes from the layout (designs with uniform layers are over-predicted,
%   designs with wide width steps under-predicted); a grade-transition term
%   is not significant. Not calibrated for RIS cables. A screening and
%   sizing model: the chosen design must be verified with WP_MECH_SURROGATE.

% coefficients (fit_jacket_surrogate.py; see the calibration note above)
JT_REF = 3e-3; BETA = 0.2;
A_M = 0.7532; B_M = 0.2151; C_M = 0.1983;   % Pm
A_B = 1.3436; B_B = 0.5489; C_B = 0.5011;   % Pm + Pb
if nargin < 8 || isempty(JT)
    JT = JT_REF*ones(size(sigma_nom));   % no JT given: no thickness correction (beta term = 1)
end
nt = n_turns(:)';
F = nt.*Iop.*B_layer(:)';
q = (cumsum(F) - 0.5*F)/W1;
sigma_acc = q/p_rs.*sigma_nom(:)';
step = max(0, nt - [nt(2:end) nt(end)])./nt;      % width step under each layer
gJ = (JT(:)'/JT_REF).^BETA;
Pm   = S_z + gJ.*(A_M*sigma_nom(:)' + B_M*sigma_acc + C_M*sigma_acc.*step);
PmPb = S_z + gJ.*(A_B*sigma_nom(:)' + B_B*sigma_acc + C_B*sigma_acc.*step);
end
