function out = cs_smoke_test(run_fem, overrides)
%CS_SMOKE_TEST Run the whole CS chain without any dialog, on the shipped test baseline.
%
%   out = CS_SMOKE_TEST() reads the machine geometry in this folder
%   (baseline_VNS_07_2026_V4_LG.csv: 6 CS modules, PF1..PF6, plasma; 7
%   scenarios), scans the CS design space with the parameters of the input
%   template and returns the results in out.DATA. With run_fem=true it also
%   verifies the first feasible design point with FEM_CS_VERIFY. overrides
%   (optional struct) replaces input parameters; the default narrows the
%   scan (30 mm minimum cable size, 5 kA current step) so it runs quickly.
%
%   Works in MATLAB (parameters from the input workbook) and in GNU Octave
%   (parameters from params_template.txt, a text copy of the template -
%   regenerate it when the template changes; test/octave_shims/ provides
%   stand-ins for TABLE/HEIGHT, which Octave lacks, and is only added to
%   the path in Octave).

if nargin < 1, run_fem = false; end
if nargin < 2, overrides = struct('min_size_CICC', 0.03, 'Iop_steps', 5000); end
here = fileparts(mfilename('fullpath'));
cs = fileparts(here);
% test/ holds Octave-only stand-ins for TABLE/HEIGHT that must never shadow MATLAB's own
addpath(fileparts(cs)); made_paths('CS');
is_octave = exist('OCTAVE_VERSION', 'builtin') ~= 0;

if is_octave
    addpath(fullfile(here, 'octave_shims'));
    p = load_params_txt(fullfile(here, 'params_template.txt'));
    d = csvread(fullfile(here, 'baseline_VNS_07_2026_V4_LG.csv'));
else
    p = read_machine_input(fullfile(cs, 'input', 'WP_CS_input_template.xlsx'));
    d = readmatrix(fullfile(here, 'baseline_VNS_07_2026_V4_LG.csv'));
end

ov = fieldnames(overrides);
for k = 1:numel(ov)
    p.(ov{k}) = overrides.(ov{k});
end

MAt = max(abs(d(:, 5:end)), [], 2);
g = compute_operating_params(p);
env = generate_combinations(p, g, MAt);
fprintf('turns %d..%d, layers %d..%d, %d combT x %d combL x %d Iop\n', ...
    env.n_turns_min, env.n_turns_max, env.min_n_layers, env.max_n_layers, ...
    size(env.combT, 1), size(env.combL, 1), numel(env.Iop_range));
DATA = scan_wp_designs(p, g, env);
fprintf('%d feasible design point(s)\n', height(DATA));

out.DATA = DATA; out.p = p; out.g = g; out.fem = [];
if run_fem && height(DATA) > 0
    geom.R = d(:, 1); geom.Z = d(:, 2); geom.dr = d(:, 3); geom.dz = d(:, 4);
    geom.MAt_signed = d(:, 5:end);
    out.fem = fem_cs_verify(DATA(1), p, g, geom);
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
