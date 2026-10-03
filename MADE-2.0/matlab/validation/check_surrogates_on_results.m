function C = check_surrogates_on_results(results_file, idx, out_csv)
%CHECK_SURROGATES_ON_RESULTS Verify the jacket and case surrogates with the 2D FE on saved scan solutions.
%
%   C = CHECK_SURROGATES_ON_RESULTS(results_file, idx) reloads the design
%   points idx (default 1:3, i.e. the smallest radial builds) of a scan
%   saved by main_WP_TF_design (results xlsx with its .mat companion), runs
%   the 2D FE (WP_MECH_SURROGATE, primary load case) on each and compares
%   the jacket and case Pm and Pm+Pb with the surrogate values of the scan
%   (columns JT_Pm, JT_PmPb, S_T_VT = case Pm, Case_PmPb). It also reports
%   whether the design is inside the calibration range (SURROGATE_VALIDITY).
%
%   C = CHECK_SURROGATES_ON_RESULTS(results_file, idx, out_csv) appends one
%   line per design to out_csv (default
%   validation/results/surrogate_checks.csv): the machine, the design and
%   FE / surrogate values. These runs are the material to extend the
%   calibration sets to other machines.
%
%   Each FE run takes about 10 minutes.

if nargin < 2 || isempty(idx), idx = 1:3; end
this_dir = fileparts(mfilename('fullpath'));
if nargin < 3 || isempty(out_csv), out_csv = fullfile(this_dir, 'results', 'surrogate_checks.csv'); end
[~, tag] = fileparts(results_file);
new_file = ~isfile(out_csv);
fid = fopen(out_csv, 'a');
if new_file
    fprintf(fid, ['results,idx,n_TF,Ri,Iop,n_layers,W1,case_x,inside,FE_valid,' ...
        'jacket_Pm_FE,jacket_Pm_sur,jacket_PmPb_FE,jacket_PmPb_sur,case_Pm_FE,case_Pm_sur,case_PmPb_FE,case_PmPb_sur\n']);
end
C = struct('idx', {}, 'FE', {}, 'sur', {}, 'ratio', {}, 'inside', {});
for i = idx(:)'
    [row, p] = load_design_point(results_file, i);
    V = surrogate_validity(row, p, false);
    r = to_struct(row);
    fe = wp_mech_surrogate(r, p, struct('classify', 1, 'verbose', false));
    P = fe.primary;
    FE = 1e-6*[max(P.layer.Pm), max(P.layer.PmPb), max([P.case_scl.Pm]), max([P.case_scl.PmPb])];
    sur = [getv(row, 'JT_Pm'), getv(row, 'JT_PmPb'), getv(row, 'S_T_VT'), getv(row, 'Case_PmPb')];
    ratio = sur./FE;
    inside = V.machine_ok && V.n_out == 0;
    C(end+1) = struct('idx', i, 'FE', FE, 'sur', sur, 'ratio', ratio, 'inside', inside); %#ok<AGROW>
    fprintf(['#%d (%s): FE valid %d | jacket Pm %.0f / %.0f, Pm+Pb %.0f / %.0f | case Pm %.0f / %.0f, ' ...
        'Pm+Pb %.0f / %.0f MPa (FE / surrogate) | surrogate/FE %s | inside calibration %d\n'], i, tag, fe.valid, ...
        FE(1), sur(1), FE(2), sur(2), FE(3), sur(3), FE(4), sur(4), mat2str(round(ratio*1000)/1000), inside);
    if ~isempty(V.lines), fprintf('%s\n', V.lines{:}); end
    fprintf(fid, '%s,%d,%d,%g,%g,%d,%g,%g,%d,%d,%g,%g,%g,%g,%g,%g,%g,%g\n', tag, i, p.n_TF, r.Ri_, r.Iop, r.n_layers, ...
        V.value(3), V.value(5), inside, fe.valid, FE(1), sur(1), FE(2), sur(2), FE(3), sur(3), FE(4), sur(4));
end
fclose(fid);
end

function v = getv(row, f)
v = NaN;
if istable(row) && any(strcmp(row.Properties.VariableNames, f)), v = row.(f); end
if isstruct(row) && isfield(row, f), v = row.(f); end
end

function s = to_struct(row)
% table row -> struct with the per-layer fields cut to n_layers
if isstruct(row), s = row; return, end
s = table2struct(row);
nl = s.n_layers;
for f = {'n_turns', 'Cond_w', 'Cond_h', 'JT', 'type_cable', 'S_Cable', 'r_cable'}
    if isfield(s, f{1}), v = s.(f{1}); s.(f{1}) = v(1:min(nl, numel(v))); end
end
end
