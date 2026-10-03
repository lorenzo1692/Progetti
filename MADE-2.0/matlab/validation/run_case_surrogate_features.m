function run_case_surrogate_features(p, csv_file)
%RUN_CASE_SURROGATE_FEATURES Inputs of the case stress surrogate fit.
%
%   RUN_CASE_SURROGATE_FEATURES(p, csv_file) evaluates the layered model
%   (WP_LAYERED_CYLINDER, default options) on the 19 reference 2D FE runs
%   (validation/results/mech_reference_rows.mat: the exact rows the FE ran,
%   RUN_MECH_REFERENCE_FE) and writes one row per run: nose thickness at the
%   centre plane t_n, WP height h_wp, WP widths, radii, hoop resultant H,
%   axial strain, layered-model nose and plate stresses, centering force.
%   validation/tools/fit_case_surrogate.py fits CASE_STRESS_SURROGATE on
%   these features and validation/results/mech_reference_fe.csv.
%   Default csv_file: validation/results/case_surrogate_features.csv.

this_dir = fileparts(mfilename('fullpath'));
if nargin < 2 || isempty(csv_file)
    csv_file = fullfile(this_dir, 'results', 'case_surrogate_features.csv');
end
L = load(fullfile(this_dir, 'results', 'mech_reference_rows.mat'));
fid = fopen(csv_file, 'w');
fprintf(fid, 'name,t_n,h_wp,W1,Wmax,nl,Ri,Rb,Rj,H,eps0,lcPm,lcPmPb,lcst,lcsz,lcplPm,F_c,nt_tot,Iop\n');
for i = 1:numel(L.ref)
    row = L.ref(i).row;
    o = wp_layered_cylinder(row, p, struct());
    Rb = o.rings(1).r1; Rj = o.rings(1).r2; Rp = o.rings(end).r1;
    H = 0;
    for k = 1:numel(o.rings), H = H + trapz(o.sample(k).r, o.sample(k).st); end
    nl = row.n_layers; W = row.n_turns(1:nl).*row.Cond_w(1:nl);
    fprintf(fid, '%s,%g,%g,%g,%g,%d,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%d,%g\n', L.ref(i).name, Rj - Rb, Rp - Rj, ...
        W(1), max(W), nl, row.Ri_, Rb, Rj, H, o.eps0, o.nose.Pm, o.nose.PmPb, o.nose.membrane(2), ...
        o.nose.membrane(3), o.plate.Pm, o.check.F_centering_per_coil, sum(row.n_turns(1:nl)), row.Iop);
end
fclose(fid);
fprintf('%d runs written to %s\n', numel(L.ref), csv_file);
end
