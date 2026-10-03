function out = case_stress_surrogate(row, p, opts)
%CASE_STRESS_SURROGATE Case primary stresses (max over the case SCLs) calibrated on the 2D FE.
%
%   out = CASE_STRESS_SURROGATE(row, p) estimates the largest primary
%   membrane stress intensity Pm and membrane + bending Pm+Pb of the TF
%   case (nose, vault diagonal, lateral walls, plasma-side plate: the SCLs
%   of WP_MECH_SURROGATE) of the design point row, for the criteria
%   Pm <= Sm_case and Pm+Pb <= 1.5 Sm_case.
%
%   Model (docs/MODELLO_A_STRATI.md, section 11):
%     Pm   = g_m(x) * Pm_nose,LC      g_m(x) = a_m + b_m ln x
%     Pm+Pb = g_b(x) * Pm_nose,LC     g_b(x) = a_b + b_b ln x
%     x = h_WP / t_nose
%   Pm_nose,LC is the membrane stress intensity of the nose of the layered
%   model (WP_LAYERED_CYLINDER, default options): it carries the loads
%   (Lorentz hoop force, vertical force) and the stiffness of nose, WP and
%   plate. The layered model is axisymmetric and misses the frame action
%   of the case (the lateral walls hang the nose to the plate): its nose
%   error grows with the ratio of the WP height h_WP (plasma-side plate to
%   nose) to the nose thickness at the centre plane t_nose, which g
%   corrects.
%
%   Calibration: 19 2D FE runs (13 designs, 6 jacket-thickness variants,
%   validation/results/mech_reference_fe.csv, x = 2.2 ... 21), fitted by
%   validation/tools/fit_case_surrogate.py. Leave-one-design-out, model/FE:
%     Pm     0.967 ... 1.040, rms 2.1 %
%     Pm+Pb  0.864 ... 1.121, rms 7.1 %
%   Validity: n_TF = 12, rectangular cables, x within the calibration
%   range; out.extrapolated flags x outside it.
%
%   The coefficients hold for the default options of WP_LAYERED_CYLINDER
%   (03/10/2026): refit (run_case_surrogate_features + fit_case_surrogate.py)
%   if they change.
%   opts (optional): lc (precomputed WP_LAYERED_CYLINDER result).
%   out: Pm, PmPb [Pa], x, t_nose, h_wp, Pm_lc, g_m, g_b, extrapolated, lc.

if nargin < 3, opts = struct(); end
C = case_surrogate_coefficients();
if isfield(opts, 'lc') && ~isempty(opts.lc)
    lc = opts.lc;
else
    lc = wp_layered_cylinder(row, p, struct());
end
Rb = lc.rings(1).r1; Rj = lc.rings(1).r2; Rp = lc.rings(end).r1;
t_nose = Rj - Rb; h_wp = Rp - Rj;
x = h_wp/max(t_nose, 1e-6);
g_m = C.a_m + C.b_m*log(x);
g_b = C.a_b + C.b_b*log(x);
out = struct('Pm', g_m*lc.nose.Pm, 'PmPb', g_b*lc.nose.Pm, 'x', x, 't_nose', t_nose, 'h_wp', h_wp, ...
    'Pm_lc', lc.nose.Pm, 'g_m', g_m, 'g_b', g_b, 'extrapolated', x < C.x_range(1) || x > C.x_range(2), 'lc', lc);
end

function C = case_surrogate_coefficients()
% fitted by validation/tools/fit_case_surrogate.py on the 19 2D FE runs
C.a_m = 0.95802119; C.b_m = 0.13122639;     % max case Pm
C.a_b = 0.99183351; C.b_b = 0.23057292;     % max case Pm+Pb
C.x_range = [2.2, 21.1];
end
