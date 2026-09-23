function plot_wp_diagnostics(row, p, fig_title)
%PLOT_WP_DIAGNOSTICS Per-layer EM/thermal companion plots for one design point.
%
%   PLOT_WP_DIAGNOSTICS(row, p) plots, layer by layer, three profiles
%   across the winding pack for one design point (row, as for
%   PLOT_WP_SECTION):
%     1. Hot-spot temperature THS against its type-specific allowable
%        (p.THS_max_LTS / p.THS_max_HTS) - where it saturates across grades.
%     2. Peak field: the smeared/linear model used during the scan
%        (B_layers) AND the exact discrete Biot-Savart field at every turn
%        (COMPUTE_DISCRETE_FIELD_PROFILE), so the ripple the smeared model
%        misses is visible directly, with the peak enhancement factor
%        printed on the plot and in the console.
%     3. Engineering current density Iop/(Cond_w*Cond_h) per layer.
%   Bars/markers are colored by cable type (LTS/HTS) so a grade change is
%   visible.

if nargin < 3 || isempty(fig_title)
    fig_title = sprintf('WP diagnostics - Iop=%.0f A, n_layers=%d', row.Iop, row.n_layers);
end

n_layers = row.n_layers;
n_turns = row.n_turns(1:n_layers);
Cond_w  = row.Cond_w(1:n_layers);
Cond_h  = row.Cond_h(1:n_layers);
THS     = row.THS(1:n_layers);
type_cable = row.type_cable(1:n_layers);
is_hts = strcmp(type_cable, 'HTS');

ths_limit = repmat(p.THS_max_LTS, 1, n_layers);
ths_limit(is_hts) = p.THS_max_HTS;

% Reconstruct the smeared field profile exactly as the "Recompute B per
% grade" section of search/scan_wp_designs.m does.
n_spire_ = zeros(1, n_layers);
n_spire_(1) = sum(n_turns);
for k = 2:n_layers
    n_spire_(k) = n_spire_(k-1) - n_turns(k);
end
B_layers = row.B_TF .* (n_spire_/n_spire_(1));

% Exact discrete-coil field (see physics/compute_discrete_field_profile.m)
field = compute_discrete_field_profile(row, p);
B_discrete_max_per_layer = zeros(1, n_layers);
for k = 1:n_layers
    B_discrete_max_per_layer(k) = max(field.B_discrete(field.layer == k));
end
[peak_ripple, peak_layer] = max(field.ripple);
fprintf(['Discrete-coil field check: peak B_discrete/B_smooth ripple = %.3f at layer %d ' ...
    '(B_smooth=%.3f T, B_discrete=%.3f T)\n'], peak_ripple, field.layer(peak_layer), ...
    field.B_smooth(peak_layer), field.B_discrete(peak_layer));

JENG_layer = row.Iop ./ (Cond_w.*Cond_h) * 1e-6; % [A/mm^2]

col_LTS = [0.20 0.45 0.85];
col_HTS = [0.85 0.35 0.15];
bar_colors = repmat(col_LTS, n_layers, 1);
bar_colors(is_hts, :) = repmat(col_HTS, sum(is_hts), 1);

figure('units', 'normalized', 'outerposition', [0.1 0.1 0.55 0.85]);
sgtitle(fig_title, 'Interpreter', 'none');

subplot(3,1,1)
b = bar(1:n_layers, THS, 'FaceColor', 'flat');
b.CData = bar_colors;
hold on
stairs((0.5:n_layers+0.5), [ths_limit, ths_limit(end)], 'r--', 'LineWidth', 1.5);
hold off
xlabel('Layer'); ylabel('THS [K]');
title('Hot-spot temperature per layer (dashed = type-specific allowable)');
grid on

subplot(3,1,2)
b = bar(1:n_layers, B_layers, 'FaceColor', 'flat');
b.CData = bar_colors;
hold on
plot(1:n_layers, B_discrete_max_per_layer, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5, ...
    'DisplayName', 'Discrete Biot-Savart (max per layer)');
hold off
xlabel('Layer'); ylabel('B [T]');
title(sprintf('Field profile: smeared model (bars) vs discrete coils (markers) - peak ripple %.2fx', peak_ripple));
grid on
legend('Location', 'eastoutside');

subplot(3,1,3)
b = bar(1:n_layers, JENG_layer, 'FaceColor', 'flat');
b.CData = bar_colors;
xlabel('Layer'); ylabel('J_{eng} [A/mm^2]');
title('Engineering current density per layer');
grid on

legend_handles = [patch(nan,nan,col_LTS), patch(nan,nan,col_HTS)];
legend(legend_handles, {'LTS','HTS'}, 'Location', 'eastoutside');
end
