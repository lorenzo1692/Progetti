function out = size_jacket_surrogate(in)
%SIZE_JACKET_SURROGATE Minimum jacket thickness per grade from the fast jacket surrogate.
%
%   out = SIZE_JACKET_SURROGATE(in) finds, for a winding pack whose layout
%   (turns per layer, cell width, cable area per layer, field per layer) is
%   fixed, the smallest jacket thickness of every grade such that the
%   primary jacket stresses of JACKET_STRESS_SURROGATE, multiplied by the
%   margin, satisfy the criteria of the 2D FE figure of merit:
%
%       margin * Pm_k    <= Sm         margin * (Pm+Pb)_k <= 1.5 Sm
%
%   in every layer k. Rectangular cables (shape_cable = 201) only.
%
%   For a given layout the surrogate stress of layer k depends only on the
%   jacket thickness of that layer: sigma_nom_k through the cell geometry
%   (same rule as SIZE_CICC_CABLE: r_SC = clamp(JT, r_SC_min, r_SC_max),
%   SC_w = Cond_w - 2 JT - 2 t_ins, SC_h = (A_cable + (4 - pi) r_SC^2)/SC_w,
%   Cond_h = SC_h + 2 JT + 2 t_ins) and the accumulated load q_k, which does
%   not depend on JT. So every layer is sized on its own, by increasing JT
%   in steps of JT_step from its starting value, and the grade takes the
%   largest JT of its layers (one conductor per grade). The axial stress
%   S_z is an input: the caller iterates with the case sizing, which
%   changes S_z (more jacket steel, different nose).
%
%   in fields (1 x n_layers vectors unless noted):
%     A_cable, Cond_w, tins, E_cbl, n_turns, B_layer, grade (grade index of
%     each layer), JT0 (starting thickness)
%     scalars: E_jckt, E_ins, r_SC_min, r_SC_max, p_rs, Iop, S_z, Sm,
%     margin, JT_step, JT_max
%
%   out: JT (per layer, constant within each grade), Pm, PmPb (surrogate
%   stresses at that JT, without the margin), ok (all layers within the
%   criteria at JT <= JT_max), crit_layer (layer with the highest Pm+Pb
%   utilization), Cond_h (cell height at that JT).

nl = numel(in.A_cable);
JT = in.JT0(:)';
W1 = in.Cond_w(1)*in.n_turns(1);
ok = false;
for it = 1:ceil((in.JT_max - min(JT))/in.JT_step) + 2
    [Pm, PmPb, Ch] = stresses(JT, in, W1);
    bad = in.margin*Pm > in.Sm | in.margin*PmPb > 1.5*in.Sm;
    if ~any(bad)
        ok = true; break
    end
    if any(JT(bad) + in.JT_step > in.JT_max + 1e-12)
        break
    end
    JT(bad) = JT(bad) + in.JT_step;
end
% one jacket thickness per grade
for g = unique(in.grade(:)')
    m = in.grade == g;
    JT(m) = max(JT(m));
end
[Pm, PmPb, Ch] = stresses(JT, in, W1);
ok = ok && all(in.margin*Pm <= in.Sm & in.margin*PmPb <= 1.5*in.Sm);
[~, crit] = max(max(Pm/in.Sm, PmPb/(1.5*in.Sm)));
out = struct('JT', JT, 'Pm', Pm, 'PmPb', PmPb, 'ok', ok, 'crit_layer', crit, 'Cond_h', Ch);
end

function [Pm, PmPb, Ch] = stresses(JT, in, W1)
r = min(max(JT, in.r_SC_min), in.r_SC_max);
tins = in.tins(:)'.*ones(size(JT));
SCw = in.Cond_w - 2*JT - 2*tins;
SCh = (in.A_cable + (4 - pi)*r.^2)./SCw;
Ch = SCh + 2*JT + 2*tins;
Ke = 2*in.E_jckt*JT./Ch + 2*tins*in.E_ins./Ch + ...
    1./(1./(in.E_cbl.*SCw./SCh) + 2./(in.E_jckt*in.Cond_w./JT) + 2./(in.E_ins*in.Cond_w./tins));
dcr_jckt = 2*in.E_jckt*JT./Ch./Ke;
r_steel = (in.Cond_w - 2*tins)./(2*JT);
sigma_nom = in.p_rs*r_steel.*dcr_jckt;
[Pm, PmPb] = jacket_stress_surrogate(sigma_nom, in.p_rs, in.n_turns, in.B_layer, in.Iop, W1, in.S_z, JT);
end
