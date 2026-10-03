function R = refine_jacket_fe(row, p, opts)
%REFINE_JACKET_FE Verify the jacket of a chosen design with the 2D FE and correct it.
%
%   R = REFINE_JACKET_FE(row, p) is the last step of the jacket sizing
%   (docs/DIMENSIONAMENTO_JACKET.md, section 7). The scan sizes the jacket
%   of every grade with the fast surrogate and a uniform jacket_margin; the
%   surrogate's error depends on the layout (median FE/surrogate 1.03,
%   90th percentile 1.17 on Pm+Pb), but its response to the jacket
%   thickness of a GIVEN layout is within ~3% of the FE. So, for the chosen
%   design only:
%
%     1. WP_MECH_SURROGATE (2D FE) gives the primary Pm and Pm+Pb of the
%        jacket per layer; if Pm <= Sm and Pm+Pb <= 1.5 Sm everywhere, stop;
%     2. otherwise the FE/surrogate ratio of each layer becomes the margin
%        of that layer (times 1 + opts.extra, default 2%), SIZE_JACKET_SURROGATE
%        re-sizes the jacket of every grade with formulas, JACKET_JT_VARIANT
%        builds the new WP (same cable areas and cell width, same nose
%        thickness: the nose moves inward with the taller WP, conservative
%        for the case, which a thicker jacket unloads);
%     3. back to 1, at most opts.max_iter FE runs (default 3).
%
%   A layer the FE finds safer than the surrogate keeps its jacket (the
%   ratio is not allowed below 1 for re-sizing, no thinning). The new WP must
%   still fit the sector (toroidal gap >= p.toroidal_gap) or the refinement
%   stops with R.ok = false and a message.
%
%   Inputs: row - design point (table row of the scan or struct with Iop,
%   n_layers, n_turns, Cond_w, Cond_h, JT, Ri_, Rk_, type_cable,
%   shape_cable, B_grade); p - input parameters; opts - max_iter, extra,
%   and any WP_MECH_SURROGATE option.
%
%   R: row (refined design, struct), ok, history (one entry per FE run: JT,
%   FE and surrogate maxima, Rk_, radial build), mech (last FE result).

if nargin < 3, opts = struct(); end
max_iter = getf(opts, 'max_iter', 3); extra = getf(opts, 'extra', 0.02);
fe_opts = rmf(opts, {'max_iter', 'extra'});
r = to_struct(row);
nl = r.n_layers;
if get_shape(r, p) ~= 201
    error('refine_jacket_fe:ris', 'Rectangular cables only (the jacket surrogate is not calibrated for RIS).');
end
Mu_0 = 4e-7*pi;
Sm = getf(p, 'Sm_jacket', p.S_amm_JT);
hist = struct('JT', {}, 'FE_Pm', {}, 'FE_PmPb', {}, 'sur_Pm', {}, 'sur_PmPb', {}, 'Rk_', {}, 'radial_build', {});
R = struct('row', r, 'ok', false, 'history', hist, 'mech', []);
for it = 1:max_iter
    mech = wp_mech_surrogate(r, p, fe_opts);
    if ~mech.valid
        fprintf(2, 'refine_jacket_fe: 2D FE not valid (%s), stop.\n', mech.checks.summary);
        R.mech = mech; R.row = r; return
    end
    FE_Pm = mech.primary.layer.Pm(1:nl); FE_PmPb = mech.primary.layer.PmPb(1:nl);
    s = surrogate_in(r, p, Mu_0, Sm);
    s.evaluate_only = true; s.margin = 1;
    ev = size_jacket_surrogate(s);
    R.history(end+1) = struct('JT', r.JT(1:nl), 'FE_Pm', max(FE_Pm), 'FE_PmPb', max(FE_PmPb), ...
        'sur_Pm', max(ev.Pm), 'sur_PmPb', max(ev.PmPb), 'Rk_', r.Rk_, 'radial_build', r.Ri_ - r.Rk_/cos(pi/p.n_TF));
    fprintf('refine_jacket_fe pass %d: JT %s mm | FE Pm %.0f, Pm+Pb %.0f MPa | surrogate Pm %.0f, Pm+Pb %.0f MPa | radial build %.1f mm\n', ...
        it, mat2str(unique(round(r.JT(1:nl)*1e4)/10)), max(FE_Pm)/1e6, max(FE_PmPb)/1e6, max(ev.Pm)/1e6, ...
        max(ev.PmPb)/1e6, 1e3*(r.Ri_ - r.Rk_/cos(pi/p.n_TF)));
    R.mech = mech; R.row = r;
    if all(FE_Pm <= Sm) && all(FE_PmPb <= 1.5*Sm)
        R.ok = true; return
    end
    if it == max_iter, break, end
    % per-layer FE/surrogate ratio as the margin of that layer, re-size
    s = rmf(s, {'evaluate_only'});
    s.margin = max(1, FE_Pm./ev.Pm)*(1 + extra);
    s.margin_PmPb = max(1, FE_PmPb./ev.PmPb)*(1 + extra);
    s.JT_step = p.JT_step; s.JT_max = getf(p, 'JT_max', 0.010);
    js = size_jacket_surrogate(s);
    if ~js.ok
        fprintf(2, 'refine_jacket_fe: no jacket up to JT_max satisfies the corrected criteria, stop.\n');
        return
    end
    v = jacket_jt_variant(r, p, js.JT - r.JT(1:nl));
    if v.min_toroidal_gap < getf(p, 'toroidal_gap', 0.015)
        fprintf(2, ['refine_jacket_fe: with the thicker jacket the WP no longer fits the sector ' ...
            '(toroidal gap %.1f mm < %.1f mm): the layout must change (fewer turns per layer), stop.\n'], ...
            1e3*v.min_toroidal_gap, 1e3*getf(p, 'toroidal_gap', 0.015));
        R.row = v; return
    end
    r = v;
end
fprintf(2, 'refine_jacket_fe: criteria not satisfied after %d FE runs.\n', max_iter);
end

function s = surrogate_in(r, p, Mu_0, Sm)
% inputs of SIZE_JACKET_SURROGATE for the design r (as the scan builds them)
nl = r.n_layers;
tg = wp_turn_geometry(r, p);
E_cbl = p.E_cbl_LTS*ones(1, nl);
for k = 1:nl
    if strcmp(r.type_cable{k}, 'HTS'), E_cbl(k) = p.E_cbl_HTS; end
end
A = tg.A_cable;
grade = cumsum([1, abs(diff(A)) > 1e-12*max(A)]);          % one grade per cable
Bl = wp_peak_field_fast(r, p);
p_rs = r.B_grade(1)^2/(2*Mu_0);
g = compute_operating_params(p);
theta = 2*pi/p.n_TF;
Ch = r.Cond_h(1:nl); Cw = r.Cond_w(1:nl); nt = r.n_turns(1:nl);
A_WP = sum(Ch.*Cw.*nt); A_JT_tot = sum(tg.A_jacket.*nt);
R_bore = r.Rk_/cos(theta/2);                               % arc bore, as SIZE_CASE_VAULT
A_CASE = r.Ri_^2*tan(theta/2) - theta/2*R_bore^2 - A_WP;
T_bf = axial_load_factor(p)*0.5*(g.k_bf*p.n_TF*(sum(nt)*r.Iop)^2*Mu_0/(2*pi));
s = struct('A_cable', A, 'Cond_w', Cw, 'tins', tg.tins, 'E_cbl', E_cbl, 'n_turns', nt, 'B_layer', Bl, ...
    'grade', grade, 'JT0', r.JT(1:nl), 'E_jckt', p.E_jckt, 'E_ins', p.E_ins, 'r_SC_min', p.r_SC_min, ...
    'r_SC_max', p.r_SC_max, 'p_rs', p_rs, 'Iop', r.Iop, 'S_z', T_bf/(A_JT_tot + A_CASE), 'Sm', Sm, ...
    'JT_step', p.JT_step, 'JT_max', getf(p, 'JT_max', 0.010));
end

function s = to_struct(row)
% table row of the scan -> struct with the per-layer fields cut to n_layers
if isstruct(row)
    s = row;
else
    s = struct();
    for v = row.Properties.VariableNames, s.(v{1}) = row.(v{1}); end
end
nl = s.n_layers;
for f = {'n_turns', 'Cond_w', 'Cond_h', 'JT', 'B_grade', 'S_Cable', 'type_cable'}
    if isfield(s, f{1}), x = s.(f{1}); s.(f{1}) = x(1:nl); end
end
end

function c = get_shape(r, p)
if isfield(r, 'shape_cable') && ~isempty(r.shape_cable), c = r.shape_cable(1); else, c = p.shape_cable; end
end

function v = getf(s, f, d)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end

function s = rmf(s, names)
for i = 1:numel(names)
    if isfield(s, names{i}), s = rmfield(s, names{i}); end
end
end
