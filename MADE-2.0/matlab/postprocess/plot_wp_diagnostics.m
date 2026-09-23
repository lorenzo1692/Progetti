function plot_wp_diagnostics(row, fig_title)
%PLOT_WP_DIAGNOSTICS Per-layer EM/thermal companion plots for one design point.
%
%   PLOT_WP_DIAGNOSTICS(row) plots, layer by layer, three profiles across
%   the winding pack for one design point (row, as for PLOT_WP_SECTION):
%     1. Hot-spot temperature THS - where it is highest/closest to
%        saturating across grades.
%     2. Peak field B_layers (reconstructed the same way
%        search/scan_wp_designs.m does, from B_TF and the turns profile).
%     3. Engineering current density Iop/(Cond_w*Cond_h) per layer.
%   Bars are colored by cable type (LTS/HTS) so a grade change is visible.
%
%   NOTE: this tool does not yet carry a hot-spot temperature allowable as
%   an input parameter, so no limit line is drawn on the THS plot - read
%   the values against your own criterion for now.

if nargin < 2 || isempty(fig_title)
    fig_title = sprintf('WP diagnostics - Iop=%.0f A, n_layers=%d', row.Iop, row.n_layers);
end

n_layers = row.n_layers;
n_turns = row.n_turns(1:n_layers);
Cond_w  = row.Cond_w(1:n_layers);
Cond_h  = row.Cond_h(1:n_layers);
THS     = row.THS(1:n_layers);
type_cable = row.type_cable(1:n_layers);

% Reconstruct the field profile exactly as the "Recompute B per grade"
% section of search/scan_wp_designs.m does.
n_spire_ = zeros(1, n_layers);
n_spire_(1) = sum(n_turns);
for k = 2:n_layers
    n_spire_(k) = n_spire_(k-1) - n_turns(k);
end
B_layers = row.B_TF .* (n_spire_/n_spire_(1));

JENG_layer = row.Iop ./ (Cond_w.*Cond_h) * 1e-6; % [A/mm^2]

col_LTS = [0.20 0.45 0.85];
col_HTS = [0.85 0.35 0.15];
bar_colors = repmat(col_LTS, n_layers, 1);
is_hts = strcmp(type_cable, 'HTS');
bar_colors(is_hts, :) = repmat(col_HTS, sum(is_hts), 1);

figure('units', 'normalized', 'outerposition', [0.1 0.1 0.55 0.85]);
sgtitle(fig_title, 'Interpreter', 'none');

subplot(3,1,1)
b = bar(1:n_layers, THS, 'FaceColor', 'flat');
b.CData = bar_colors;
xlabel('Layer'); ylabel('THS [K]');
title('Hot-spot temperature per layer');
grid on

subplot(3,1,2)
b = bar(1:n_layers, B_layers, 'FaceColor', 'flat');
b.CData = bar_colors;
xlabel('Layer'); ylabel('B [T]');
title('Peak field profile across the WP');
grid on

subplot(3,1,3)
b = bar(1:n_layers, JENG_layer, 'FaceColor', 'flat');
b.CData = bar_colors;
xlabel('Layer'); ylabel('J_{eng} [A/mm^2]');
title('Engineering current density per layer');
grid on

legend_handles = [patch(nan,nan,col_LTS), patch(nan,nan,col_HTS)];
legend(legend_handles, {'LTS','HTS'}, 'Location', 'eastoutside');
end
