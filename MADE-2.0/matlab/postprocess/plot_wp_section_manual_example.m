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

% Field on the WP (BSUM peak): needed by plot_wp_diagnostics and to size
% the cable of each grade below.
row.B_TF = 13.489;               % [T]
row.Tau_discharge = 20;          % [s] discharge time constant

% Conductor grades: first layer of each grade. Each grade's cable is sized
% by cicc() at the smeared field of its first layer (as the scan does);
% set row.N_Sc / row.N_Cu / row.type_cable by hand instead to test a
% specific cable.
grade_start = [1 7];
n_spire_ = sum(row.n_turns) - [0 cumsum(row.n_turns(2:end))];
B_layers = row.B_TF * n_spire_/n_spire_(1);
row.N_Sc = zeros(1, row.n_layers); row.N_Cu = zeros(1, row.n_layers);
row.THS = zeros(1, row.n_layers);  row.B_grade = zeros(1, row.n_layers);
for k = grade_start
    [tc, N_Cu, N_Sc, ~, ~, ~, ~, THS] = cicc(B_layers(k), row.Iop, row.Tau_discharge, ...
        p.WP_SC_type, p.THS_max_LTS, p.THS_max_HTS);
    row.type_cable(k:end) = tc; row.N_Sc(k:end) = N_Sc; row.N_Cu(k:end) = N_Cu;
    row.THS(k:end) = THS;       row.B_grade(k:end) = B_layers(k);
end

%% Plot
plot_wp_section(row, p, 'Manual what-if configuration');
plot_wp_section_bfield(row, p, 'Manual what-if configuration');
plot_wp_diagnostics(row, p, 'Manual what-if configuration');
plot_hotspot_transient(row, p, 'Manual what-if configuration - hot spot');

% 2D FE mechanical surrogate (true equilibrium, rounded turns, contacts)
mech = wp_mech_surrogate(row, p);
plot_wp_mech_surrogate(mech, p, 'Manual what-if configuration');
