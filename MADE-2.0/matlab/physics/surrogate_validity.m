function V = surrogate_validity(DATA, p, verbose)
%SURROGATE_VALIDITY Are the scan solutions inside the calibration range of the FE surrogates?
%
%   V = SURROGATE_VALIDITY(DATA, p) checks the machine (p) and every design
%   point of DATA (table of SCAN_WP_DESIGNS, or struct array) against the
%   range of the 2D FE runs on which the jacket surrogate
%   (JACKET_STRESS_SURROGATE, scf_model = 1) and the case surrogate
%   (CASE_STRESS_SURROGATE, case_model = 1) were calibrated:
%     n_TF = 12, rectangular cables, Ri_ 1.18-1.26 m, Iop 21-66 kA,
%     toroidal WP width at the plasma side W1 0.30-0.54 m, 11-28 layers,
%     case x = h_WP/t_nose 2.2-21.
%   Outside these ranges the surrogates extrapolate: the chosen design must
%   be verified with the 2D FE (validation/check_surrogates_on_results.m),
%   and the new FE runs can be added to the calibration sets.
%
%   V.machine_ok (n_TF, shape), V.out (per design: logical per quantity),
%   V.n_out (designs outside the range for at least one quantity),
%   V.names, V.range, V.lines (report). verbose (default true) prints the report.

if nargin < 3, verbose = true; end
R = struct('Ri', [1.18 1.26], 'Iop', [21e3 66e3], 'W1', [0.30 0.545], 'n_layers', [11 28], 'case_x', [2.2 21.1]);
names = fieldnames(R)';
if isstruct(DATA), n = numel(DATA); else, n = height(DATA); end
val = nan(n, numel(names));
for i = 1:n
    if isstruct(DATA), r = DATA(i); else, r = DATA(i, :); end
    val(i, 1) = r.Ri_; val(i, 2) = r.Iop;
    nt = r.n_turns; cw = r.Cond_w;
    val(i, 3) = nt(1)*cw(1); val(i, 4) = r.n_layers;
    if has(r, 'case_x'), val(i, 5) = r.case_x; end
end
out = false(n, numel(names));
for j = 1:numel(names)
    lim = R.(names{j});
    out(:, j) = ~isnan(val(:, j)) & (val(:, j) < lim(1) | val(:, j) > lim(2));
end
shape = 201; if isfield(p, 'shape_cable'), shape = p.shape_cable; end
V = struct('machine_ok', p.n_TF == 12 && shape == 201, 'names', {names}, 'range', R, 'value', val, ...
    'out', out, 'n_out', sum(any(out, 2)), 'lines', {{}});
L = {};
if p.n_TF ~= 12
    L{end+1} = sprintf('  machine: n_TF = %d (calibrated on 12): every solution must be verified with the 2D FE', p.n_TF);
end
if shape ~= 201
    L{end+1} = sprintf('  machine: shape_cable = %d (calibrated on 201, Rect): the scan uses the analytic formulas', shape);
end
for j = 1:numel(names)
    if any(out(:, j))
        lim = R.(names{j}); v = val(:, j);
        L{end+1} = sprintf('  %-9s %4d of %d solutions outside [%g, %g] (min %g, max %g)', names{j}, ...
            sum(out(:, j)), n, lim(1), lim(2), min(v), max(v)); %#ok<AGROW>
    end
end
V.lines = L;
if verbose
    if isempty(L)
        fprintf('Surrogates: all %d solutions inside the calibration range of the jacket and case surrogates.\n', n);
    else
        fprintf(['Surrogates: %d of %d solutions outside the calibration range of the jacket/case surrogates ' ...
            '(extrapolation; verify the chosen design with the 2D FE):\n'], V.n_out, n);
        fprintf('%s\n', L{:});
    end
end
end

function t = has(r, f)
if isstruct(r), t = isfield(r, f); else, t = any(strcmp(r.Properties.VariableNames, f)); end
end
