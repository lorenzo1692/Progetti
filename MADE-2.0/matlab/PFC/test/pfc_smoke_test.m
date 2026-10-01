function out = pfc_smoke_test(n_PF, WP_h, run_fem, overrides)
%PFC_SMOKE_TEST Run the whole PFC chain without any dialog, on one coil at one WP_h.
%
%   out = PFC_SMOKE_TEST(n_PF, WP_h, run_fem, overrides) reads the machine geometry
%   shipped in this folder (baseline_VNS_07_2026_V4_LG.csv, 13 rows: 6 CS,
%   PF1..PF6, plasma; 7 scenarios), computes the scenario forces and the
%   coupling matrix, scans the PF coil n_PF (default 6) at winding-pack
%   height WP_h (default 0.30 m), and, with run_fem=true (default), verifies
%   the first feasible design point with FEM_PFC_VERIFY. Returns the
%   results table (struct array in Octave) in out.DATA and the FEM result
%   in out.fem. overrides (optional struct) replaces input parameters, e.g.
%   struct('fz_source', 1). Takes about a minute.
%
%   Works in MATLAB (parameters from the input workbook) and in GNU Octave
%   (parameters from params_template.txt, a text copy of the template -
%   regenerate it when the template changes; test/octave_shims/ provides
%   stand-ins for TABLE/HEIGHT, which Octave lacks, and is only added to
%   the path in Octave).

if nargin < 1, n_PF = 6; end
if nargin < 2, WP_h = 0.30; end
if nargin < 3, run_fem = true; end
if nargin < 4, overrides = struct(); end

here = fileparts(mfilename('fullpath'));
pfc = fileparts(here);
addpath(genpath(pfc)); addpath(fileparts(pfc));
is_octave = exist('OCTAVE_VERSION', 'builtin') ~= 0;

if is_octave
    addpath(fullfile(here, 'octave_shims'));
    p = load_params_txt(fullfile(here, 'params_template.txt'));
    d = csvread(fullfile(here, 'baseline_VNS_07_2026_V4_LG.csv'));
else
    p = read_machine_input(fullfile(pfc, 'input', 'WP_PFC_input_template.xlsx'));
    d = readmatrix(fullfile(here, 'baseline_VNS_07_2026_V4_LG.csv'));
end

ov = fieldnames(overrides);
for k = 1:numel(ov)
    p.(ov{k}) = overrides.(ov{k});
end

geom.R = d(:, 1); geom.Z = d(:, 2); geom.dr = d(:, 3); geom.dz = d(:, 4);
geom.MAt_signed = d(:, 5:end); geom.MAt_scenario = abs(geom.MAt_signed);
geom.MAt = max(geom.MAt_scenario, [], 2);

Fz = compute_scenario_forces(geom);
FZ_max = max(abs(Fz), [], 2)*1e6;
L = compute_coupling_matrix(geom);
g = compute_operating_params(p);

row = p.cs_modules + n_PF;
L_factor = sum(L(:, row))/L(row, row);
env = generate_combinations(p, g, WP_h);
DATA = scan_wp_designs(p, g, env, WP_h, geom.R(row), geom.MAt(row), FZ_max(row), L_factor, n_PF);
fprintf('PF%d, WP_h = %.3f m: %d feasible design point(s)\n', n_PF, WP_h, height(DATA));

out.DATA = DATA; out.fem = [];
if run_fem && height(DATA) > 0
    out.fem = fem_pfc_verify(DATA(1, :), p, geom);
end
end

function p = load_params_txt(fn)
fid = fopen(fn); p = struct();
while true
    l = fgetl(fid);
    if ~ischar(l), break; end
    t = strsplit(strtrim(l), ' ');
    p.(t{1}) = str2double(t{2});
end
fclose(fid);
end
