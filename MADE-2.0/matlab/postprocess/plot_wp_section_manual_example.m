%% plot_wp_section_manual_example.m
% Template for a manual, hand-edited "what-if" WP configuration - no scan,
% no Excel results file. Copy this script, edit the numbers below, and run
% it to preview a cross-section (e.g. to try a specific FEM-derived
% geometry directly against physics/size_case_vault.m's assumptions).
%
% Machine-level parameters (materials, allowables, n_TF, ...) still come
% from a real input file, since those aren't something you'd normally
% hand-edit per plot; only the winding-pack geometry below is manual.

clearvars; clc
this_dir = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(this_dir, '..')));

p = read_machine_input(fullfile(this_dir, '..', 'input', 'WP_TF_input_template.xlsx'));

%% Manual WP configuration - edit freely
row = struct();
row.Iop = 62300;                                   % [A] operating current
row.n_layers = 11;                                 % number of radial layers (rows)
row.n_turns  = [10 10 10 10 10 10 8 8 8 8 8];       % turns per layer
row.Cond_w   = repmat(0.046, 1, 11);                % [m] cell width (toroidal), per layer
row.Cond_h   = [repmat(0.0326,1,6) repmat(0.0228,1,5)]; % [m] cell height (radial), per layer
row.JT       = repmat(0.0035, 1, 11);               % [m] jacket thickness, per layer
row.type_cable = repmat({'LTS'}, 1, 11);            % 'LTS' or 'HTS', per layer

row.Ri_ = 1.18;   % [m] Case outer (plasma-side) radius
row.Rj_ = 0.8354; % [m] WP inner boundary radius
row.Rk_ = 0.70;   % [m] Case nose tip radius

% Only needed for plot_wp_diagnostics (hot-spot/field/current-density
% profile) - fill these in from a cicc() run or a saved DATA row if you
% want that plot too; plot_wp_section alone does not need them.
row.B_TF = 13.489;               % [T] peak field on the WP (BSUM)
row.THS  = [65 66 67 68 69 70 72 73 74 75 76]; % [K] hot-spot temperature per layer

%% Plot
plot_wp_section(row, p, 'Manual what-if configuration');
plot_wp_diagnostics(row, p, 'Manual what-if configuration');
