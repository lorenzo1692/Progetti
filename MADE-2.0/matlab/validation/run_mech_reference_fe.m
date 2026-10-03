function run_mech_reference_fe(p, out_dir, worker, n_workers, csv_file)
%RUN_MECH_REFERENCE_FE 2D FE reference results (jacket and case) of the reference designs.
%
%   RUN_MECH_REFERENCE_FE(p, out_dir, worker, n_workers) runs
%   WP_MECH_SURROGATE (2D FE, validated against ANSYS) on the reference
%   designs of validation/results/mech_reference_designs.mat (the 13 base
%   designs of the jacket surrogate calibration) and on jacket-thickness
%   variants of three of them, and saves one file out_dir/fe_<name>.mat per
%   run with, for the primary (Lorentz + axial) and total (+ cool-down)
%   load cases: jacket Pm and Pm+Pb per layer, case SCL Pm and Pm+Pb, the
%   axial strain and the validity checks. Runs already on disk are skipped
%   (restartable); jobs are split among n_workers parallel sessions.
%
%   RUN_MECH_REFERENCE_FE(p, out_dir, 0, 0, csv_file) collects the runs into
%   csv_file (validation/results/mech_reference_fe.csv): one row per run
%   and SCL / layer, long format (name, kind, index, label, Pm_P, PmPb_P,
%   Pm_T, PmPb_T), kind = 'layer' or 'scl'.
%
%   These results are the reference of the layered-cylinder model
%   (physics/wp_layered_cylinder.m, docs/MODELLO_A_STRATI.md): the case SCLs
%   of the primary load case were not stored by the jacket calibration runs.

if nargin < 3 || isempty(worker), worker = 1; end
if nargin < 4 || isempty(n_workers), n_workers = 1; end
this_dir = fileparts(mfilename('fullpath'));
if ~exist(out_dir, 'dir'), [~] = mkdir(out_dir); end
B = load(fullfile(this_dir, 'results', 'mech_reference_designs.mat'));
jobs = struct('name', {}, 'row', {}, 'dJT', {});
for b = 1:numel(B.base)
    jobs(end+1) = struct('name', B.base(b).name, 'row', B.base(b).row, 'dJT', 0); %#ok<AGROW>
end
for b = 1:numel(B.base)          % jacket-thickness variants: jacket/case coupling
    if any(strcmp(B.base(b).name, {'bench', 'd7', 'd10'}))
        for dj = [0.5 1.0]
            jobs(end+1) = struct('name', sprintf('%s@jt+%.1f', B.base(b).name, dj), 'row', B.base(b).row, ...
                'dJT', dj*1e-3); %#ok<AGROW>
        end
    end
end
if worker == 0
    collect(jobs, out_dir, csv_file);
    return
end
p.shape_cable = 201;
for j = worker:n_workers:numel(jobs)
    f = fullfile(out_dir, ['fe_' safe(jobs(j).name) '.mat']);
    if exist(f, 'file'), continue, end
    row = jobs(j).row;
    if jobs(j).dJT > 0, row = jacket_jt_variant(row, p, jobs(j).dJT); end
    fprintf('[%s] start\n', jobs(j).name);
    try
        out = wp_mech_surrogate(row, p, struct('classify', 1));
        r = struct('name', jobs(j).name, 'row', row, 'valid', out.valid, 'checks', out.checks.summary, ...
            'eps_z', out.eps_z, 'layer_Pm_P', out.primary.layer.Pm, 'layer_PmPb_P', out.primary.layer.PmPb, ...
            'layer_Pm_T', out.layer.Pm, 'layer_PmPb_T', out.layer.PmPb, ...
            'scl_name', {{out.case_scl.name}}, 'scl_Pm_P', [out.primary.case_scl.Pm], ...
            'scl_PmPb_P', [out.primary.case_scl.PmPb], 'scl_Pm_T', [out.case_scl.Pm], 'scl_PmPb_T', [out.case_scl.PmPb]);
        save('-v7', f, 'r');
        fprintf('[%s] done: %s\n', jobs(j).name, out.checks.summary);
    catch err
        fprintf('[%s] FAILED: %s\n', jobs(j).name, err.message);
    end
end
end

function collect(jobs, out_dir, csv_file)
fid = fopen(csv_file, 'w');
fprintf(fid, 'name,kind,index,label,Pm_P,PmPb_P,Pm_T,PmPb_T,valid\n');
n = 0;
for j = 1:numel(jobs)
    f = fullfile(out_dir, ['fe_' safe(jobs(j).name) '.mat']);
    if ~exist(f, 'file'), continue, end
    S = load(f); r = S.r; n = n + 1;
    for k = 1:numel(r.layer_Pm_P)
        fprintf(fid, '%s,layer,%d,L%d,%.6e,%.6e,%.6e,%.6e,%d\n', r.name, k, k, r.layer_Pm_P(k), r.layer_PmPb_P(k), ...
            r.layer_Pm_T(k), r.layer_PmPb_T(k), r.valid);
    end
    for k = 1:numel(r.scl_Pm_P)
        fprintf(fid, '%s,scl,%d,%s,%.6e,%.6e,%.6e,%.6e,%d\n', r.name, k, r.scl_name{k}, r.scl_Pm_P(k), r.scl_PmPb_P(k), ...
            r.scl_Pm_T(k), r.scl_PmPb_T(k), r.valid);
    end
end
fclose(fid);
fprintf('%d run(s) written to %s\n', n, csv_file);
end

function s = safe(name)
s = regexprep(name, '[^A-Za-z0-9_+.-]', '_');
end
