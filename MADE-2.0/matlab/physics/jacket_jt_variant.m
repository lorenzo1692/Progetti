function v = jacket_jt_variant(row, p, dJT)
%JACKET_JT_VARIANT Same winding pack with a different jacket thickness.
%
%   v = JACKET_JT_VARIANT(row, p, dJT) returns the design point row with the
%   jacket thickness of every layer increased by dJT [m] (scalar, or one
%   value per layer), keeping what the scan keeps when it sizes the jacket
%   of a given layout:
%     - the superconducting cable area of each layer (A_cable from
%       WP_TURN_GEOMETRY, i.e. the cable the conductor sizing produced),
%     - the cell width Cond_w (turns per layer x Cond_w = WP width),
%     - the case nose thickness DTF = Rj_ - Rk_ (the nose moves radially
%       inward with the taller WP: Rk_ is recomputed).
%   The cable height and the cell height follow from the area:
%     SC_w = Cond_w - 2 JT - 2 t_ins, r_SC = clamp(JT, r_SC_min, r_SC_max),
%     SC_h = (A_cable + (4 - pi) r_SC^2)/SC_w, Cond_h = SC_h + 2 JT + 2 t_ins.
%   v.min_toroidal_gap is the smallest lateral distance between the
%   grounded WP and the case flank over the layers (the quantity the scan
%   requires to be >= p.toroidal_gap): the variant does not shrink the turns
%   as the scan would, so a thick jacket can push the deep layers out of
%   the sector (gap < 0).
%   Rectangular cables only. Used to build the jacket-thickness variants of
%   the calibration set of JACKET_STRESS_SURROGATE
%   (validation/run_jacket_jt_calibration.m).

nl = row.n_layers;
tg = wp_turn_geometry(row, p);
if tg.is_round
    error('jacket_jt_variant:ris', 'Rectangular cables only.');
end
dJT = dJT(:)'.*ones(1, nl);
tins = tg.tins;
rmin = getf(p, 'r_SC_min', 2e-3); rmax = getf(p, 'r_SC_max', 6e-3);
Cw = row.Cond_w(1:nl);
A = tg.A_cable;
JT = row.JT(1:nl) + dJT;
r = min(max(JT, rmin), rmax);
SCw = Cw - 2*JT - 2*tins;
SCh = (A + (4 - pi)*r.^2)./SCw;
Ch = SCh + 2*JT + 2*tins;
ins = getf(p, 'INS_grades', 5e-4); dps = getf(p, 'dr_plasma_side', 0.02); git = getf(p, 'GoundIns', 5e-3);
Rj_old = row.Ri_ - sum(row.Cond_h(1:nl)) - (nl-1)*ins - dps - 2*git;
Rj_new = row.Ri_ - sum(Ch) - (nl-1)*ins - dps - 2*git;
v = row;
v.JT(1:nl) = JT;
v.Cond_h(1:nl) = Ch;
v.Rk_ = row.Rk_ - (Rj_old - Rj_new);
if isfield(v, 'S_Cable'), v.S_Cable(1:nl) = A; end
Re = row.Ri_ - dps - git; gap = zeros(1, nl);
for k = 1:nl
    Ri = Re - Ch(k);
    gap(k) = (2*Ri*tan(pi/getf(p, 'n_TF', 12)) - (Cw(k)*row.n_turns(k) + 2*git))/2;
    Re = Ri - ins;
end
v.min_toroidal_gap = min(gap);
end

function x = getf(s, f, d)
if isfield(s, f) && ~isempty(s.(f)), x = s.(f); else, x = d; end
end
