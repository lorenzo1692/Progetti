%% validate_mech_surrogate_2026.m
% 2D FE mechanical surrogate (physics/wp_mech_surrogate.m) vs the two ANSYS
% 2D generalized-plane-strain runs of 2026:
%   A. TF_FEM_benchmark_2026 (21-Sep): 11 layers (6x10 + 5x8 turns),
%      Cond_w 46 mm, Cond_h 32.6/22.8 mm, JT 3.5 mm, r_SC 4 mm, 62.3 kA
%   B. design point 7 (EM_2D007, 24-Sep): 13 layers x 8 turns, 3 grades,
%      Cond_w 43.0 mm, Cond_h 34.2/26.1/24.2 mm, JT 3.0 mm, 64.628 kA
%
% FEM reference values below were extracted from the TFBM_*.txt exports
% (nodal-averaged, material-restricted stresses), linearized on exactly the
% same sections the surrogate uses (jacket walls and fillets, case SCLs),
% with the scripts in validation/tools (see the README there).
% See TF_FEM_benchmark_2026_findings.md, "2D FE mechanical surrogate".
%
% Acceptance (review C03): every quantity must fall inside the band in
% TOL below and the surrogate result must pass its own validity checks;
% the script ends with an error if anything fails. The axial force T_bf is
% prescribed from each FEM run (to compare like with like); the value the
% tool would compute from the machine input is printed next to it.

clearvars; clc
this_dir = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(this_dir, '..')));
p = read_machine_input(fullfile(this_dir, '..', 'input', 'WP_TF_input_template.xlsx'));

TOL.eps_z = 0.02;                 % |relative error| on the axial strain
TOL.peak  = [0.85 1.10];          % surrogate/FEM, jacket peak per layer
TOL.Pm    = [0.90 1.10];          % surrogate/FEM, max linearized Pm per layer
TOL.PmPb  = [0.85 1.10];          % surrogate/FEM, max linearized Pm+Pb per layer
TOL.scl   = [0.90 1.05; 0.90 1.05; 0.90 1.05; 0.95 1.05; 0.95 1.05; 0.95 1.05; 0.95 1.20];  % case SCL Pm
failures = {};

cases = struct([]);
% --- A. benchmark ------------------------------------------------------
c.name = 'TF_FEM_benchmark_2026';
c.row = struct('Iop', 62300, 'n_layers', 11, 'n_turns', [10 10 10 10 10 10 8 8 8 8 8], ...
    'Cond_w', repmat(0.046,1,11), 'Cond_h', [repmat(0.0326,1,6) repmat(0.0228,1,5)], ...
    'JT', repmat(0.0035,1,11), 'Ri_', 1.18, 'Rk_', 0.70);
c.row.type_cable = repmat({'LTS'}, 1, 11);
c.opts = struct('T_bf', 32675062.45, 'r_SC', 0.004);
c.eps_z = -1.28176e-3;
c.peak = [613 690 742 783 839 980 695 749 798 846 876];                    % max SINT per layer [MPa]
c.Pm   = [487 472 485 503 526 581 505 516 530 543 549];                    % max linearized Pm per layer
c.PmPb = [632 640 648 693 739 863 634 677 717 754 778];                    % max linearized Pm+Pb per layer
c.scl  = [697 736; 689 700; 670 779; 418 420; 411 413; 389 400; 492 508];  % case SCLs [Pm PmPb]
cases = [cases, c];
% --- B. design 7 ---------------------------------------------------------
c.name = 'design 7 (EM_2D007)';
c.row = struct('Iop', 64628, 'n_layers', 13, 'n_turns', 8*ones(1,13), ...
    'Cond_w', repmat(0.04300426499905924,1,13), ...
    'Cond_h', [repmat(0.034167368089179834,1,4) repmat(0.026104457386404455,1,3) repmat(0.024226622311504957,1,6)], ...
    'JT', repmat(0.003,1,13), 'Ri_', 1.259425287356322, 'Rk_', 0.7300827089713595);
c.row.type_cable = repmat({'LTS'}, 1, 13);
c.opts = struct('T_bf', 36.266e6);
c.eps_z = -1.57384e-3;
c.peak = [674 864 932 977 879 929 958 939 965 984 1004 1025 1124];
c.Pm   = [474 459 481 508 503 524 540 543 552 560 568 573 581];
c.PmPb = [641 718 791 840 772 815 841 827 844 856 868 889 963];
c.scl  = [639 745; 639 657; 647 789; 363 370; 347 354; 335 339; 505 532];
cases = [cases, c];

for i = 1:numel(cases)
    c = cases(i);
    fprintf('\n===== %s =====\n', c.name);
    out = wp_mech_surrogate(c.row, p, c.opts);
    g = compute_operating_params(p);
    T_tool = 0.5*g.k_bf*p.n_TF*(sum(c.row.n_turns)*c.row.Iop)^2*(4e-7*pi)/(2*pi);
    fprintf('axial force T_bf: prescribed from the FEM %.3f MN (tool formula with this input: %.3f MN)\n', ...
        c.opts.T_bf/1e6, T_tool/1e6);
    fprintf('validity checks: %s\n', out.checks.summary);
    if ~out.valid, failures{end+1} = sprintf('%s: validity checks failed (%s)', c.name, out.checks.summary); end %#ok<SAGROW>
    fprintf('axial strain eps_z: FEM %.4e, surrogate %.4e (%+.1f%%)\n', c.eps_z, out.eps_z, 100*(out.eps_z/c.eps_z-1));
    if abs(out.eps_z/c.eps_z - 1) > TOL.eps_z, failures{end+1} = sprintf('%s: eps_z', c.name); end %#ok<SAGROW>
    fprintf('%-6s %16s %16s %16s   [MPa, FEM/surrogate]\n', 'layer', 'jacket peak', 'max Pm', 'max Pm+Pb');
    for k = 1:numel(c.peak)
        fprintf('L%-5d %7.0f/%-7.0f  %7.0f/%-7.0f  %7.0f/%-7.0f\n', k, c.peak(k), out.layer.peak(k)/1e6, ...
            c.Pm(k), out.layer.Pm(k)/1e6, c.PmPb(k), out.layer.PmPb(k)/1e6);
    end
    r = [out.layer.peak/1e6./c.peak; out.layer.Pm/1e6./c.Pm; out.layer.PmPb/1e6./c.PmPb];
    fprintf('ratio surrogate/FEM (mean, min, max): peak %.3f %.3f %.3f | Pm %.3f %.3f %.3f | Pm+Pb %.3f %.3f %.3f\n', ...
        mean(r(1,:)), min(r(1,:)), max(r(1,:)), mean(r(2,:)), min(r(2,:)), max(r(2,:)), mean(r(3,:)), min(r(3,:)), max(r(3,:)));
    fprintf('global jacket peak: FEM %.0f, surrogate %.0f MPa\n', max(c.peak), out.fom.jacket_peak/1e6);
    lims = [TOL.peak; TOL.Pm; TOL.PmPb]; nm = {'peak','Pm','Pm+Pb'};
    for q = 1:3
        bad = find(r(q,:) < lims(q,1) | r(q,:) > lims(q,2));
        if ~isempty(bad)
            failures{end+1} = sprintf('%s: %s outside [%.2f %.2f] at layer(s) %s', c.name, nm{q}, ...
                lims(q,1), lims(q,2), sprintf('%d ', bad)); %#ok<SAGROW>
        end
    end
    fprintf('case SCLs %-20s %14s %14s\n', '', 'Pm FEM/sur', 'Pm+Pb FEM/sur');
    for q = 1:numel(out.case_scl)
        fprintf('  %-28s %6.0f/%-6.0f %6.0f/%-6.0f\n', out.case_scl(q).name, c.scl(q,1), out.case_scl(q).Pm/1e6, ...
            c.scl(q,2), out.case_scl(q).PmPb/1e6);
        rr = out.case_scl(q).Pm/1e6/c.scl(q,1);
        if rr < TOL.scl(q,1) || rr > TOL.scl(q,2)
            failures{end+1} = sprintf('%s: case SCL "%s" Pm ratio %.2f outside [%.2f %.2f]', c.name, ...
                out.case_scl(q).name, rr, TOL.scl(q,1), TOL.scl(q,2)); %#ok<SAGROW>
        end
    end
    plot_wp_mech_surrogate(out, p, c.name);
end

fprintf('\n===== validation result =====\n');
if isempty(failures)
    fprintf('PASSED: all quantities inside the acceptance bands, all validity checks passed.\n');
else
    fprintf('%s\n', failures{:});
    error('validate_mech_surrogate_2026:failed', '%d validation check(s) failed.', numel(failures));
end
