%% validate_TF_FEM_design7_2026.m
% Analytical model vs ANSYS FEM for scan design point 7 (EM_2D007 run,
% 24-Sep-2026: 13 layers x 8 turns, 3 grades, Iop = 64.628 kA).
%
% Compares, at the geometry the FEM was built on:
%   1. peak field on the conductor: smeared B_TF used by the scan vs the
%      discrete Biot-Savart model (physics/compute_discrete_field_profile.m)
%   2. Jacket (per layer) and Case stress: analytical formulas
%      (validation/forward_eval_wp_stress.m) vs FEM nodal SINT.
% Findings: validation/TF_FEM_benchmark_2026_findings.md, "Design 7".

clearvars; clc
this_dir = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(this_dir, '..')));

%% Design point 7 (row of the scan results file, WP_TF_input_template_design_7.xlsx)
row = struct();
row.Iop = 64628;
row.n_layers = 13;
row.n_turns = 8*ones(1,13);
row.Cond_w = repmat(0.04300426499905924, 1, 13);
row.Cond_h = [repmat(0.034167368089179834,1,4) repmat(0.026104457386404455,1,3) repmat(0.024226622311504957,1,6)];
row.JT = [repmat(0.0031,1,4) repmat(0.0030,1,9)];
row.Ri_ = 1.259425287356322;
row.Rk_ = 0.7300827089713595;
row.B_TF = 13.482710702701707;            % smeared peak used by the scan
row.lateral_w = 0.15374619886744612;
row.S_T_JT = 651.1; row.S_T_VT = 666.2;   % [MPa] scan prediction
B_grade_scan = [13.483 9.334 6.223];      % field each grade was sized at
grade_first = [1 5 8];

p = struct('n_TF',12,'dr_plasma_side',0.02,'GoundIns',0.005,'INS_grades',0.0005, ...
    'turn_insulation_nominal',0.001,'Increm',1);

%% FEM reference (TFBM_*.txt nodal-averaged, EM_2D002.png)
FEM_B_peak = 14.482;                      % [T] BSUM max on conductor
FEM_JT_SINT = [673.7 863.7 932.0 977.4 879.0 928.8 957.6 939.4 965.3 984.4 1003.6 1025.3 1124.4]; % [MPa] max per layer
FEM_JT_peak_node = 169508;                % L13, outermost turn, cable fillet corner next to nose/side wall
FEM_Case_SINT = 800.5;                    % [MPa] node 184060, nose
T_bf = 36.266e6;                          % [N] GSBDATA axial force per leg (CHECK_TF / T_BF/2)
% NOTE: the FEM mesh uses JT = 3.0 mm in grade 1 (cable 35.0 x 26.2 mm), the
% scan row has 3.1 mm; the FEM geometry is used below.

%% 1. Field
row_fem = row; row_fem.JT(:) = 0.003;
field = compute_discrete_field_profile(row_fem, p);
B_disc = max(field.B_discrete);
fprintf('Peak field on conductor: smeared %.3f T | discrete %.3f T | FEM %.3f T\n', row.B_TF, B_disc, FEM_B_peak);
fprintf('Field each grade was sized at vs discrete peak on that grade:\n');
for g = 1:3
    k = grade_first(g);
    fprintf('  grade %d (from L%d): sized at %.2f T, discrete %.2f T (%+.0f%%)\n', g, k, ...
        B_grade_scan(g), max(field.B_discrete(field.layer == k)), ...
        100*(max(field.B_discrete(field.layer == k))/B_grade_scan(g)-1));
end

%% 2. Stress
geo = struct('n_turns',row.n_turns,'Cond_h',row.Cond_h,'Cond_w',row.Cond_w,'JT',row_fem.JT, ...
    'tins',0.001*ones(1,13),'Ri_',row.Ri_,'Rk_',row.Rk_,'lateral_w',row.lateral_w,'n_TF',12, ...
    'dr_plasma_side',0.02,'r_SC_min',0.002,'r_SC_max',0.006);
geo.Rj_ = row.Ri_ - sum(row.Cond_h) - 12*p.INS_grades - p.dr_plasma_side - 2*p.GoundIns; % FEM WP_INTR
mat = struct('E_jckt',205e9,'E_case',205e9,'E_ins',12e9,'E_cbl',10e9);
is_tr = false(1,13); is_tr([4 5 7 8]) = true;

cases = {row.B_TF, 3.15, 'B smeared, transition SCF 3.15 (= scan)'; ...
         B_disc,   3.15, 'B discrete, transition SCF 3.15'; ...
         B_disc,   1.00, 'B discrete, no SCF'};
fprintf('\n%-42s %9s %9s\n', 'Analytical at FEM geometry', 'Jacket', 'Case');
for c = 1:size(cases,1)
    o = forward_eval_wp_stress(geo, mat, cases{c,1}, T_bf, is_tr, cases{c,2});
    fprintf('%-42s %7.0f M %7.0f M\n', cases{c,3}, o.S_T_JT_max/1e6, o.S_T_VT/1e6);
    if c == 3, o_noscf = o; end
end
fprintf('%-42s %7.0f M %7.0f M\n', 'FEM (SINT)', max(FEM_JT_SINT), FEM_Case_SINT);
fprintf('%-42s %7.0f M %7.0f M\n', 'scan row', row.S_T_JT, row.S_T_VT);

fprintf('\nJacket per layer [MPa]: FEM SINT vs analytical (B discrete, no SCF)\n');
for k = 1:13
    fprintf('  L%-2d FEM %6.0f   analytical %5.0f   ratio %.2f\n', k, FEM_JT_SINT(k), ...
        o_noscf.S_T_JT(k)/1e6, FEM_JT_SINT(k)/(o_noscf.S_T_JT(k)/1e6));
end
