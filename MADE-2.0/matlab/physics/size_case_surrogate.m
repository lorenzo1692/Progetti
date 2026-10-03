function out = size_case_surrogate(row, p, opts)
%SIZE_CASE_SURROGATE Size the case nose with the case stress surrogate.
%
%   out = SIZE_CASE_SURROGATE(row, p, opts) finds the smallest nose thickness
%   at the centre plane, DTF = DTF0 + k*DTF_step (k = 0, 1, ...), such that
%   CASE_STRESS_SURROGATE gives
%     margin   x Pm    <= Sm_case
%     margin_b x Pm+Pb <= 1.5 Sm_case.
%   The case stresses decrease monotonically with DTF, so k is found by
%   bisection (about 10 evaluations of the layered model, < 0.1 s).
%   Bore = arc (as SIZE_CASE_VAULT): R_bore = Rj_ - DTF, Rk_ = R_bore cos(pi/n_TF).
%
%   row: Iop, n_layers, n_turns, Cond_w, Cond_h, JT, Ri_, type_cable,
%   shape_cable (Rk_ is set here). opts: DTF0, DTF_step, DTF_max (default
%   0.6*Rj_), Sm (Sm_case), margin (case_margin, 1.05), margin_b
%   (case_margin_PmPb, 1.10).
%   out: ok, DTF, R_bore, Rk_, Rj_, Pm, PmPb [Pa], x, extrapolated, evaluations.

if nargin < 3, opts = struct(); end
nl = row.n_layers;
ins = p.INS_grades; git = p.GoundIns; dps = p.dr_plasma_side;
Rj = row.Ri_ - sum(row.Cond_h(1:nl)) - (nl - 1)*ins - dps - 2*git;
cb = cos(pi/p.n_TF);
DTF0 = getf(opts, 'DTF0', p.DTF_initial); step = getf(opts, 'DTF_step', p.DTF_step);
DTF_max = getf(opts, 'DTF_max', 0.6*Rj);
Sm = getf(opts, 'Sm', getf(p, 'Sm_case', p.S_amm_VT));
m = getf(opts, 'margin', getf(p, 'case_margin', 1.05));
mb = getf(opts, 'margin_b', getf(p, 'case_margin_PmPb', 1.10));
n_eval = 0;
check = @(k) check_k(row, p, (Rj - (DTF0 + k*step))*cb, m, mb, Sm);
kmax = floor((DTF_max - DTF0)/step);
[ok_hi, cs_hi] = check(kmax); n_eval = n_eval + 1;
if ~ok_hi
    out = struct('ok', false, 'DTF', DTF0 + kmax*step, 'R_bore', Rj - DTF0 - kmax*step, ...
        'Rk_', (Rj - DTF0 - kmax*step)*cb, 'Rj_', Rj, 'Pm', cs_hi.Pm, 'PmPb', cs_hi.PmPb, 'x', cs_hi.x, ...
        'extrapolated', cs_hi.extrapolated, 'evaluations', n_eval);
    return
end
[ok_lo, cs_lo] = check(0); n_eval = n_eval + 1;
if ok_lo
    k = 0; cs = cs_lo;
else
    lo = 0; hi = kmax; cs = cs_hi;          % check(lo) false, check(hi) true
    while hi - lo > 1
        mid = floor((lo + hi)/2);
        [okm, csm] = check(mid); n_eval = n_eval + 1;
        if okm, hi = mid; cs = csm; else, lo = mid; end
    end
    k = hi;
end
DTF = DTF0 + k*step;
out = struct('ok', true, 'DTF', DTF, 'R_bore', Rj - DTF, 'Rk_', (Rj - DTF)*cb, 'Rj_', Rj, ...
    'Pm', cs.Pm, 'PmPb', cs.PmPb, 'x', cs.x, 'extrapolated', cs.extrapolated, 'evaluations', n_eval);
end

function [ok, cs] = check_k(row, p, Rk, m, mb, Sm)
row.Rk_ = Rk;
cs = case_stress_surrogate(row, p);
ok = m*cs.Pm <= Sm && mb*cs.PmPb <= 1.5*Sm;
end

function v = getf(s, f, d)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
