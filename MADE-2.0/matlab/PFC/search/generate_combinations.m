function env = generate_combinations(p, g, WP_h)
%GENERATE_COMBINATIONS Feasible turns/layers ranges and candidate layouts for one PFC coil at one WP_h.
%
%   env = GENERATE_COMBINATIONS(p, g, WP_h) bounds the number of turns per
%   layer (from the min/max CICC size and the coil's net height WP_h -
%   2*grins_h) and the number of radial layers (from a crude
%   infinite-solenoid field window [min_B, max_B] at the extreme
%   Iop/turns), then builds every descending turns-per-grade and
%   layers-per-grade combination (combT, combL) and the operating-current
%   scan range (Iop_range) used by SCAN_WP_DESIGNS.
%
%   Relocated from the PF legacy archive (PF_opt_VNS.m, the "compute
%   bounds" and "pancake winding parameters" blocks, run once per PF coil
%   per WP_h value in the legacy driver's own double loop - see
%   MAIN_WP_PFC_DESIGN.m, which reproduces that same double loop, per the
%   29/09/2026 decision to keep PF_opt_VNS.m's all-6-coils-plus-WP_h-sweep
%   structure, rather than TF/CS's one-machine-per-run convention (see
%   manuale PFC).

available_height = WP_h - 2*p.grins_h;
env.n_turns_max = min(ceil(available_height/p.min_size_CICC), 50);
env.n_turns_min = ceil(available_height/p.max_size_CICC);

coil_height_total = WP_h*p.n_moduli + (p.n_moduli-1)*p.spacer;
env.n_layers_max = min( ...
    ceil(p.max_B*coil_height_total/(g.Mu_0*p.Iop_min*env.n_turns_min*p.n_moduli)), 40);
env.n_layers_min = ceil(p.min_B*coil_height_total/(g.Mu_0*p.Iop_max*env.n_turns_max*p.n_moduli));

Turns0 = 2*round(env.n_turns_min/2) : 2 : 2*floor(env.n_turns_max/2); % even number of turns
Layer0 = 2:round(env.n_layers_max/p.n_grades);

if isempty(Turns0) || isempty(Layer0)
    error('generate_combinations:empty_search_space', ...
        ['No feasible turns/layers range for WP_h=%.4g m (check Machine ' ...
         'geometry and Numerical settings in the input file).'], WP_h);
end

env.combT = unique(sort(nchoosek(Turns0, p.n_grades), 2, 'descend'), 'rows');
Layer_all = repmat(Layer0, 1, p.n_grades);
env.combL = unique(sort(nchoosek(Layer_all, p.n_grades), 2, 'descend'), 'rows');

env.Iop_range = p.Iop_min:p.Iop_step:p.Iop_max;
end
