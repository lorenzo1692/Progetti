function plot_wp_section(row, p, fig_title)
%PLOT_WP_SECTION Draw the 2D WP+Case cross-section for one design point.
%
%   PLOT_WP_SECTION(row, p) draws the local nose cross-section in the same
%   convention as the FEM benchmark (X = toroidal, Y = radial, see
%   benchmark_case.json "coordinates"): Case trapezoid, ground insulation,
%   every WP layer with its inter-layer insulation gap, and every turn
%   cell resolved into jacket steel / turn insulation / cable, at the
%   exact geometry this solver computed for that design point. Use it to
%   sanity-check a design visually before EXPORT_ANSYS_INPUT, or to
%   preview a manually edited "what-if" configuration.
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
%   here is this solver's own assumption, not a FEM cross-check. The
%   wedge-side insulation (FEM's CASE_INS) has no equivalent in this
%   solver's formulas and is not drawn.

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

col_case   = [0.85 0.85 0.85];
col_ground = [0.55 0.80 0.55];   % ground insulation (green)
col_layer_ins = [0.55 0.80 0.55]; % inter-layer insulation, same family as ground
col_jacket = [0.55 0.55 0.60];   % steel
col_tins   = [0.95 0.90 0.60];   % turn insulation (pale yellow)
col_LTS    = [0.20 0.45 0.85];
col_HTS    = [0.85 0.35 0.15];

% Case outline (straight-sided trapezoid - this solver's own assumption)
CASE_w_top = 2*row.Ri_*tan(theta_TF/2);
CASE_w_bot = 2*row.Rk_*tan(theta_TF/2);
case_x = [-CASE_w_top/2, CASE_w_top/2, CASE_w_bot/2, -CASE_w_bot/2];
case_y = [row.Ri_, row.Ri_, row.Rk_, row.Rk_];
h_case = patch(case_x, case_y, col_case, 'EdgeColor', [0.3 0.3 0.3], 'DisplayName', 'Case');

% Ground insulation: plasma-side wall to WP outer edge (top), and WP inner
% edge to the case nose interior (bottom). Widths follow the local case
% envelope at that radius, matching env.CASE_w's formula.
plasma_wall_y = row.Ri_ - p.dr_plasma_side;
w_top_out = 2*plasma_wall_y*tan(theta_TF/2);
w_top_in  = 2*Re(1)*tan(theta_TF/2);
patch([-w_top_out/2, w_top_out/2, w_top_in/2, -w_top_in/2], ...
      [plasma_wall_y, plasma_wall_y, Re(1), Re(1)], col_ground, 'EdgeColor', 'none');

w_bot_out = 2*Ri(n_layers)*tan(theta_TF/2);
w_bot_in  = 2*(Ri(n_layers)-p.GoundIns)*tan(theta_TF/2);
h_ground = patch([-w_bot_out/2, w_bot_out/2, w_bot_in/2, -w_bot_in/2], ...
      [Ri(n_layers), Ri(n_layers), Ri(n_layers)-p.GoundIns, Ri(n_layers)-p.GoundIns], ...
      col_ground, 'EdgeColor', 'none', 'DisplayName', 'Ground insulation');

% Inter-layer insulation gaps
h_layer_ins = [];
for k = 1:n_layers-1
    w1 = 2*Ri(k)*tan(theta_TF/2);
    w2 = 2*Re(k+1)*tan(theta_TF/2);
    r = patch([-w1/2, w1/2, w2/2, -w2/2], [Ri(k), Ri(k), Re(k+1), Re(k+1)], ...
        col_layer_ins, 'EdgeColor', 'none');
    if isempty(h_layer_ins), h_layer_ins = r; end
end

h_jacket_leg = []; h_tins_leg = []; h_lts = []; h_hts = [];

for k = 1:n_layers
    y0 = Ri(k);
    h  = Cond_h(k);
    w_total = Cond_w(k)*n_turns(k);
    x0 = -w_total/2;
    if strcmp(type_cable{k}, 'HTS')
        col = col_HTS;
    else
        col = col_LTS;
    end
    sc_w = Cond_w(k) - 2*JT(k) - 2*tins;
    sc_h = h - 2*JT(k) - 2*tins;
    for t = 1:n_turns(k)
        xc = x0 + (t-1)*Cond_w(k);
        % Jacket steel (full cell)
        rj = rectangle('Position', [xc, y0, Cond_w(k), h], ...
            'FaceColor', col_jacket, 'EdgeColor', [0.2 0.2 0.2]);
        if isempty(h_jacket_leg), h_jacket_leg = rj; end
        % Turn insulation ring (inset by JT)
        if Cond_w(k)-2*JT(k) > 0 && h-2*JT(k) > 0
            rt = rectangle('Position', [xc+JT(k), y0+JT(k), Cond_w(k)-2*JT(k), h-2*JT(k)], ...
                'FaceColor', col_tins, 'EdgeColor', 'none');
            if isempty(h_tins_leg), h_tins_leg = rt; end
        end
        % Cable (inset by JT+tins)
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
if ~isempty(h_ground)
    legend_handles(end+1) = h_ground;
    legend_labels{end+1} = 'Ground/layer insulation';
end
if ~isempty(h_jacket_leg)
    legend_handles(end+1) = patch(nan, nan, col_jacket, 'EdgeColor', [0.2 0.2 0.2]);
    legend_labels{end+1} = 'Jacket (steel)';
end
if ~isempty(h_tins_leg)
    legend_handles(end+1) = patch(nan, nan, col_tins, 'EdgeColor', 'none');
    legend_labels{end+1} = 'Turn insulation';
end
if ~isempty(h_lts)
    legend_handles(end+1) = patch(nan, nan, col_LTS, 'EdgeColor', 'none');
    legend_labels{end+1} = 'LTS cable';
end
if ~isempty(h_hts)
    legend_handles(end+1) = patch(nan, nan, col_HTS, 'EdgeColor', 'none');
    legend_labels{end+1} = 'HTS cable';
end
legend(legend_handles, legend_labels, 'Location', 'eastoutside');
grid on
hold off
end
