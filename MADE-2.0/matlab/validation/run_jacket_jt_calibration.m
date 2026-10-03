function run_jacket_jt_calibration(p, out_dir, worker, n_workers, csv_file)
%RUN_JACKET_JT_CALIBRATION 2D FE runs of jacket-thickness variants for the jacket surrogate.
%
%   RUN_JACKET_JT_CALIBRATION(p, out_dir, worker, n_workers) runs
%   WP_MECH_SURROGATE (through JACKET_SURROGATE_FE_DATA) on the
%   jacket-thickness variants of six base designs and saves one file
%   out_dir/fe_<name>.mat per run. Runs already on disk are skipped, so the
%   function can be restarted after an interruption; jobs are split among
%   n_workers parallel sessions (worker = 1..n_workers).
%
%   RUN_JACKET_JT_CALIBRATION(p, out_dir, 0, 0, csv_file) only collects the
%   runs found in out_dir and appends them to csv_file (the calibration
%   data validation/results/jacket_surrogate_calibration.csv), one row per
%   layer, in the format of the existing rows. Rows whose name is already
%   in the file are not written again.
%
%   Why: in the original calibration set every design has a single jacket
%   thickness, so the fit cannot separate the effect of JT from the effect
%   of the layout; on design 10 the surrogate's response to JT came out
%   ~1.5x stronger than the FE (docs/DIMENSIONAMENTO_JACKET.md).
%
%   Base designs (validation/results/jacket_calibration_base_rows.mat):
%   TF_FEM_benchmark_2026 (bench), design 7 (d7), design 10 (d10) and three
%   designs of the calibrated-field scan already in the set: s1_1 (14
%   layers of 8 turns, no width steps), s1_13 (25 layers, 12 -> 10 -> 8
%   turns), s2_4 (24 layers, 8 -> 6 -> 4 turns). Variants: JT + 0.5 mm and
%   JT + 1.0 mm in every layer (JACKET_JT_VARIANT: same cable area, same
%   cell width, same nose thickness); for d10 also the graded variant
%   +0 / +0.5 / +1.0 mm on its three grades. Variant names are
%   '<base>@jt+0.5', '<base>@jt+1.0', 'd10@jt+0/0.5/1.0': the part before
%   '@' identifies the base design (leave-one-design-out groups).
%   The axial force T_bf is the tool's (JACKET_SURROGATE_FE_DATA), as in the
%   existing rows.

if nargin < 3 || isempty(worker), worker = 1; end
if nargin < 4 || isempty(n_workers), n_workers = 1; end
this_dir = fileparts(mfilename('fullpath'));
if ~exist(out_dir, 'dir'), [~] = mkdir(out_dir); end   % no error if a parallel worker made it first
B = load(fullfile(this_dir, 'results', 'jacket_calibration_base_rows.mat'));
jobs = job_list(B.base);

if worker == 0
    collect(jobs, out_dir, csv_file);
    return
end
p.shape_cable = 201;
for j = worker:n_workers:numel(jobs)
    f = fullfile(out_dir, ['fe_' safe(jobs(j).name) '.mat']);
    if exist(f, 'file'), continue, end
    v = jacket_jt_variant(jobs(j).row, p, jobs(j).dJT);
    fprintf('[%s] JT %s mm, Rk_ %.4f\n', jobs(j).name, mat2str(round(unique(v.JT(1:v.n_layers))*1e4)/10), v.Rk_);
    try
        d = jacket_surrogate_fe_data(v, p, jobs(j).name, struct('verbose', false));
        d.JT = v.JT(1:v.n_layers); d.Rk_ = v.Rk_; %#ok<STRNU>
        save('-v7', f, 'd');
    catch err
        fprintf('[%s] FAILED: %s\n', jobs(j).name, err.message);
    end
end
end

function jobs = job_list(base)
jobs = struct('name', {}, 'row', {}, 'dJT', {});
for b = 1:numel(base)
    for dj = [0.5 1.0]
        jobs(end+1) = struct('name', sprintf('%s@jt+%.1f', base(b).name, dj), 'row', base(b).row, 'dJT', dj*1e-3); %#ok<AGROW>
    end
    if strcmp(base(b).name, 'd10')
        r = base(b).row; nl = r.n_layers;
        % grades of design 10: layers 1-5, 6-8, 9-13
        dJT = [zeros(1,5) 0.5*ones(1,3) 1.0*ones(1,nl-8)]*1e-3;
        jobs(end+1) = struct('name', 'd10@jt+0/0.5/1.0', 'row', r, 'dJT', dJT); %#ok<AGROW>
    end
end
end

function collect(jobs, out_dir, csv_file)
txt = '';
if exist(csv_file, 'file'), txt = fileread(csv_file); end
fid = fopen(csv_file, 'a');
n = 0;
for j = 1:numel(jobs)
    f = fullfile(out_dir, ['fe_' safe(jobs(j).name) '.mat']);
    if ~exist(f, 'file') || ~isempty(strfind(txt, [jobs(j).name ','])), continue, end
    S = load(f); d = S.d;
    for k = 1:d.n_layers
        fprintf(fid, '%s,%d,%d,%.5f,%.1f,%.5f,%.5f,%d,%.5f,%.5f,%.5f,%.5f,%.6e,%.6e,%.5f,%.6e,%.6e,%.5f,%d,%d,%.4f,%.6e,%.6e\n', ...
            d.name, k, d.n_layers, d.W, d.Iop, d.param(k), d.A(k), d.T(k), d.JT_Ch(k), d.Ch_Cw(k), d.r_steel(k), ...
            d.dcr_jckt(k), d.sigma_nom(k), d.S_z, d.dcr_WP_rad, d.PmPb_P(k), d.PmPb_T(k), d.SCF_FE(k), d.valid, ...
            d.n_turns(k), d.Bl(k), d.Pm_P(k), d.Pm_T(k));
    end
    n = n + 1;
end
fclose(fid);
fprintf('%d run(s) appended to %s\n', n, csv_file);
end

function s = safe(name)
s = regexprep(name, '[^A-Za-z0-9_+.-]', '_');
end
