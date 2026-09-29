function env = generate_combinations(p, g, MAt)
%GENERATE_COMBINATIONS Feasible turns/layers ranges and candidate layouts.
%
%   env = GENERATE_COMBINATIONS(p, g, MAt) bounds the number of turns per
%   layer (from the min/max CICC size) and the number of radial layers
%   (from the winding-pack radius, the required Ampere-turns MAt, and the
%   target field window), then builds every descending turns-per-grade and
%   layers-per-grade combination (COMBT, COMBL) and the operating-current
%   scan range (Iop_range) used by SCAN_WP_DESIGNS.
%
%   MAt is the per-location maximum Ampere-turns from READ_SCENARIO_CURRENTS;
%   only max(MAt(1:6)) is used, matching the legacy driver (see manuale
%   CS, "prime 6 scenari").
%
%   Relocated from the CS legacy archive (CS_opt_VNS.m, "Parameter
%   Initialization" + "Loop through combinations" setup). Dropped
%   n_turns_grade/n_layers_grade: computed in the original but never used
%   downstream (same class of cleanup as the TF port, see manuale CS).

env = struct();

env.n_turns_max = ceil((g.WP_h - 2*g.grins_h) / p.min_size_CICC);
env.n_turns_min = ceil((g.WP_h - 2*g.grins_h) / p.max_size_CICC);

n1_min = ceil((g.Re - g.Re/2) / p.max_size_CICC);
n1_max = ceil((g.Re - g.Re/2) / p.min_size_CICC);

MAt6 = max(MAt(1:min(6, numel(MAt))));
n2_min = ceil(MAt6 / (env.n_turns_max * p.Iop_max));
n2_max = ceil(MAt6 / (env.n_turns_min * p.Iop_min));

n3_min = ceil(p.min_B * g.h_stack / (g.Mu_0 * p.Iop_max * env.n_turns_max * p.n_moduli));
n3_max = ceil(p.max_B * g.h_stack / (g.Mu_0 * p.Iop_min * env.n_turns_min * p.n_moduli));

env.min_n_layers = max([n1_min, n2_min, n3_min]);
env.max_n_layers = min([n1_max, n2_max, n3_max]);

if p.use_pancake
    env.Turns0 = 2*(floor(env.n_turns_min/2):floor(env.n_turns_max/2));  % even turns
    env.Layer0 = 2:floor(env.max_n_layers/p.n_grades);                   % minimum 2 layers
else
    env.Turns0 = env.n_turns_min:env.n_turns_max;
    env.Layer0 = 2:2:floor(env.max_n_layers/p.n_grades);                 % even layers
end

if isempty(env.Turns0) || isempty(env.Layer0)
    error('generate_combinations:empty_search_space', ...
        ['No feasible turns/layers range for the given input parameters ', ...
         '(check Machine geometry and Numerical settings in the input file).']);
end

env.combT = unique(sort(nchoosek(env.Turns0, p.n_grades), 2, 'descend'), 'rows');
env.combL = unique(sort(nchoosek(env.Layer0, p.n_grades), 2, 'descend'), 'rows');

env.Iop_range = p.Iop_min:p.Iop_steps:p.Iop_max;
end
