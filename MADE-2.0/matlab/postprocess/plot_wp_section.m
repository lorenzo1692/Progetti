function plot_wp_section(row, p, fig_title)
%PLOT_WP_SECTION Draw the 2D WP+Case cross-section for one design point.
%
%   PLOT_WP_SECTION(row, p) draws the local nose cross-section in the same
%   convention as the FEM benchmark (X = toroidal, Y = radial, see
%   benchmark_case.json "coordinates"): Case trapezoid, every WP layer,
%   individual turn cells (jacket + cable inset), at the exact geometry
%   this solver computed for that design point. Use it to sanity-check a
%   design visually before EXPORT_ANSYS_INPUT, or to preview a manually
%   edited "what-if" configuration.
%
%   row - one row of a DATA table (from SCAN_WP_DESIGNS, in memory or
%         reloaded with LOAD_DESIGN_POINT from a saved *_results_*.xlsx),
%         OR a plain struct with the same field names for a manual,
%         hand-edited configuration (see plot_wp_section_manual_example.m).
%         Required fields: Iop, n_layers, n_turns, Cond_w, Cond_h, JT,
%         type_cable, Ri_, Rj_, Rk_.
%   p   - machine parameter struct (from READ_MACHINE_INPUT). Required
%         fields: n_TF, dr_plasma_side, GoundIns, INS_grades,
%         turn_insulation_nominal, Increm.
%
%   NOTE: this solver assumes a straight-sided trapezoidal case nose; the
%   real FEM nose is a curved arc (see
%   validation/TF_FEM_benchmark_2026_findings.md) - the Case outline drawn
%   here is this solver's own assumption, not a FEM cross-check.

if nargin < 3 || isempty(fig_title)
    fig_title = sprintf('WP section - Iop=%.0f A, n_layers=%d, Rk=%.4f m', row.Iop, row.n_layers, row.Rk_);
end

n_layers = row.n_layers;
n_turns  = row.n_turns(1:n_layers);
Cond_w   = row.Cond_w(1:n_layers);
Cond_h   = row.Cond_h(1:n_layers);
JT       = row.JT(1:n_layers);
type_cable = row.type_cable(1:n_layers);

theta_TF = 2*pi/p.n_TF;
tins = p.turn_insulation_nominal*p.Increm;

% Recompute per-layer radial positions exactly as search/scan_wp_designs.m
Re = zeros(1, n_layers+1);
Re(1) = row.Ri_ - p.dr_plasma_side - p.GoundIns;
Ri = zeros(1, n_layers);
for k = 1:n_layers
    Ri(k) = Re(k) - Cond_h(k);
    Re(k+1) = Ri(k) - p.INS_grades;
end

figure('units', 'normalized', 'outerposition', [0.1 0.1 0.7 0.8]);
hold on
axis equal

% Case outline (straight-sided trapezoid - this solver's own assumption)
CASE_w_top = 2*row.Ri_*tan(theta_TF/2);
CASE_w_bot = 2*row.Rk_*tan(theta_TF/2);
case_x = [-CASE_w_top/2, CASE_w_top/2, CASE_w_bot/2, -CASE_w_bot/2];
case_y = [row.Ri_, row.Ri_, row.Rk_, row.Rk_];
h_case = patch(case_x, case_y, [0.85 0.85 0.85], 'EdgeColor', [0.3 0.3 0.3], 'DisplayName', 'Case');

% Ground insulation boundary (dashed reference line at the WP's outer edge)
plot([-1 1]*(Re(1)*tan(theta_TF/2)), [Re(1) Re(1)], 'k--', 'HandleVisibility', 'off');

color_LTS = [0.20 0.45 0.85];
color_HTS = [0.85 0.35 0.15];
h_jacket = []; h_lts = []; h_hts = [];

for k = 1:n_layers
    y0 = Ri(k);
    h  = Cond_h(k);
    w_total = Cond_w(k)*n_turns(k);
    x0 = -w_total/2;
    if strcmp(type_cable{k}, 'HTS')
        col = color_HTS;
    else
        col = color_LTS;
    end
    sc_w = Cond_w(k) - 2*JT(k) - 2*tins;
    sc_h = h - 2*JT(k) - 2*tins;
    for t = 1:n_turns(k)
        xc = x0 + (t-1)*Cond_w(k);
        r = rectangle('Position', [xc, y0, Cond_w(k), h], ...
            'FaceColor', [0.55 0.55 0.6], 'EdgeColor', [0.2 0.2 0.2]);
        if isempty(h_jacket), h_jacket = r; end
        if sc_w > 0 && sc_h > 0
            rc = rectangle('Position', [xc+JT(k)+tins, y0+JT(k)+tins, sc_w, sc_h], ...
                'FaceColor', col, 'EdgeColor', 'none');
            if strcmp(type_cable{k}, 'HTS')
                if isempty(h_hts), h_hts = rc; end
            else
                if isempty(h_lts), h_lts = rc; end
            end
        end
    end
    text(x0 - 0.01, y0 + h/2, sprintf('L%d', k), 'FontSize', 8, 'HorizontalAlignment', 'right');
end

xlabel('Toroidal x [m]');
ylabel('Radial y [m]');
title(fig_title, 'Interpreter', 'none');

legend_handles = h_case;
legend_labels = {'Case'};
if ~isempty(h_jacket)
    legend_handles(end+1) = patch(nan, nan, [0.55 0.55 0.6], 'EdgeColor', [0.2 0.2 0.2]);
    legend_labels{end+1} = 'Jacket';
end
if ~isempty(h_lts)
    legend_handles(end+1) = patch(nan, nan, color_LTS, 'EdgeColor', 'none');
    legend_labels{end+1} = 'LTS cable';
end
if ~isempty(h_hts)
    legend_handles(end+1) = patch(nan, nan, color_HTS, 'EdgeColor', 'none');
    legend_labels{end+1} = 'HTS cable';
end
legend(legend_handles, legend_labels, 'Location', 'eastoutside');
grid on
hold off
end
