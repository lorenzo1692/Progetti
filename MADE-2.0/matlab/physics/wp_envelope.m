function env = wp_envelope(p, g)
%WP_ENVELOPE Bounding envelope for the winding pack and feasible turns/layers ranges.
%
%   env = WP_ENVELOPE(p, g) uses the raw parameters p and the derived
%   operating parameters g (from COMPUTE_OPERATING_PARAMS) to compute the
%   TF case width, the feasible winding-pack toroidal width range, and the
%   turns/layers ranges that bound the combinatorial search performed by
%   SCAN_WP_DESIGNS.

env = struct();
env.CASE_w = 2*g.R_TF_Innerleg*tan(g.theta_TF/2);            % Case width
env.lateral_w_min = env.CASE_w*p.case_wedge_frac_min;        % Case thickness, wedge side
env.lateral_w_max = env.CASE_w*p.case_wedge_frac_max;

env.WP_w_max = 2*(g.R_TF_Innerleg - p.dr_plasma_side - p.GoundIns)*tan(g.theta_TF/2) - 2*env.lateral_w_min;
env.WP_w_min = 2*(g.R_TF_Innerleg - p.dr_plasma_side - p.GoundIns)*tan(g.theta_TF/2) - 2*env.lateral_w_max;

env.Ntlmax = ceil(g.NI/p.Iop_min);
env.Ntlmin = ceil(g.NI/p.Iop_max);

env.n_turns_max = ceil((env.WP_w_max - p.GoundIns*2)/p.min_size_CICC); % Maximum feasible turns per grade
env.n_turns_min = ceil((env.WP_w_min - p.GoundIns*2)/p.max_size_CICC); % Minimum feasible turns per grade

% See fix #2 (original script history): the bound on the number of radial
% layers is only a heuristic estimate of the WP's maximum radial build
% (WP_radial_build_max_est), used solely to keep the combinatorial search
% finite. Physically invalid combinations are discarded later by the
% geometric and Tresca checks in SCAN_WP_DESIGNS.
env.max_n_layers = min( ...
    ceil((p.WP_radial_build_max_est - p.GoundIns*2)/p.min_size_CICC), ...
    floor(env.Ntlmax/env.n_turns_min));
env.min_n_layers = max( ...
    ceil((p.WP_radial_build_max_est - p.GoundIns*2)/p.max_size_CICC), ...
    ceil(env.Ntlmin/env.n_turns_max));

env.n_spire_fsbl = env.Ntlmin:env.Ntlmax;

n_turns_min_even = env.n_turns_min;
if mod(n_turns_min_even, 2) == 1
    n_turns_min_even = n_turns_min_even - 1;
end
env.turns_comb = n_turns_min_even:2:env.n_turns_max;
env.layers_comb = env.min_n_layers:min(env.max_n_layers, p.maxdim);

if isempty(env.turns_comb) || isempty(env.layers_comb)
    error('wp_envelope:empty_search_space', ...
        ['No feasible turns/layers combination for the given input parameters ', ...
         '(check Machine geometry and Numerical settings in the input file).']);
end
end
