function plot_wp_section_bfield(row, p, fig_title)
%PLOT_WP_SECTION_BFIELD WP+Case cross-section, cable colored by local field.
%
%   PLOT_WP_SECTION_BFIELD(row, p) reproduces the field-colored section
%   plot from the original Plot_WP_TF.m (Case trapezoid in gray, turn
%   insulation/jacket rings, cable region colored by a jet colormap), but
%   colors each individual turn by the peak discrete Biot-Savart field on
%   that turn's cable, self-field included (COMPUTE_DISCRETE_FIELD_PROFILE), instead of the per-grade
%   smeared B_layers value the original used - every turn gets its own
%   color, not just its layer's.
%
%   row, p - see PLOT_WP_SECTION.

if nargin < 3 || isempty(fig_title)
    fig_title = sprintf('WP section - B field, Iop=%.0f A, n_layers=%d', row.Iop, row.n_layers);
end

n_layers = row.n_layers;
n_turns  = row.n_turns(1:n_layers);
Cond_w   = row.Cond_w(1:n_layers);
Cond_h   = row.Cond_h(1:n_layers);
JT       = row.JT(1:n_layers);

theta_TF = 2*pi/p.n_TF;
tins = p.turn_insulation_nominal*p.Increm;

Re = zeros(1, n_layers+1);
Re(1) = row.Ri_ - p.dr_plasma_side - p.GoundIns;
Ri = zeros(1, n_layers);
for k = 1:n_layers
    Ri(k) = Re(k) - Cond_h(k);
    Re(k+1) = Ri(k) - p.INS_grades;
end

field = compute_discrete_field_profile(row, p);
B_min = min(field.B_discrete);
B_max = max(field.B_discrete);
if B_max == B_min, B_max = B_min + eps; end
cmap = jet(256);

figure('units', 'normalized', 'outerposition', [0.1 0.1 0.7 0.8]);
hold on
axis equal

col_case   = [0.7 0.7 0.7]; % matches the original Plot_WP_TF.m Case fill
col_ground = [0.55 0.80 0.55];
col_jacket = [0.5 0.5 0.5];
col_tins   = [0.0588 1.0 1.0];

CASE_w_top = 2*row.Ri_*tan(theta_TF/2);
CASE_w_bot = 2*row.Rk_*tan(theta_TF/2);
patch([-CASE_w_top/2, CASE_w_top/2, CASE_w_bot/2, -CASE_w_bot/2], ...
      [row.Ri_, row.Ri_, row.Rk_, row.Rk_], col_case, 'EdgeColor', [0.3 0.3 0.3]);

plasma_wall_y = row.Ri_ - p.dr_plasma_side;
w_top_out = 2*plasma_wall_y*tan(theta_TF/2);
w_top_in  = 2*Re(1)*tan(theta_TF/2);
patch([-w_top_out/2, w_top_out/2, w_top_in/2, -w_top_in/2], ...
      [plasma_wall_y, plasma_wall_y, Re(1), Re(1)], col_ground, 'EdgeColor', 'none');
w_bot_out = 2*Ri(n_layers)*tan(theta_TF/2);
w_bot_in  = 2*(Ri(n_layers)-p.GoundIns)*tan(theta_TF/2);
patch([-w_bot_out/2, w_bot_out/2, w_bot_in/2, -w_bot_in/2], ...
      [Ri(n_layers), Ri(n_layers), Ri(n_layers)-p.GoundIns, Ri(n_layers)-p.GoundIns], ...
      col_ground, 'EdgeColor', 'none');
for k = 1:n_layers-1
    w1 = 2*Ri(k)*tan(theta_TF/2);
    w2 = 2*Re(k+1)*tan(theta_TF/2);
    patch([-w1/2, w1/2, w2/2, -w2/2], [Ri(k), Ri(k), Re(k+1), Re(k+1)], col_ground, 'EdgeColor', 'none');
end

turn_i = 0;
for k = 1:n_layers
    y0 = Ri(k);
    h  = Cond_h(k);
    w_total = Cond_w(k)*n_turns(k);
    x0 = -w_total/2;
    sc_w = Cond_w(k) - 2*JT(k) - 2*tins;
    sc_h = h - 2*JT(k) - 2*tins;
    for t = 1:n_turns(k)
        turn_i = turn_i + 1;
        xc = x0 + (t-1)*Cond_w(k);
        rectangle('Position', [xc, y0, Cond_w(k), h], 'FaceColor', col_tins, 'EdgeColor', [0.3 0.3 0.3]);
        if Cond_w(k)-2*tins > 0 && h-2*tins > 0
            rectangle('Position', [xc+tins, y0+tins, Cond_w(k)-2*tins, h-2*tins], ...
                'FaceColor', col_jacket, 'EdgeColor', 'none');
        end
        if sc_w > 0 && sc_h > 0
            ci = 1 + round((field.B_discrete(turn_i)-B_min)/(B_max-B_min)*255);
            rectangle('Position', [xc+JT(k)+tins, y0+JT(k)+tins, sc_w, sc_h], ...
                'FaceColor', cmap(ci,:), 'EdgeColor', 'none');
        end
    end
    text(x0 - 0.01, y0 + h/2, sprintf('L%d', k), 'FontSize', 8, 'HorizontalAlignment', 'right');
end

colormap(cmap);
set(gca, 'CLim', [B_min B_max]);
c = colorbar;
ylabel(c, 'Peak B on cable [T] (discrete Biot-Savart + self-field)');

xlabel('Toroidal x [m]');
ylabel('Radial y [m]');
title(fig_title, 'Interpreter', 'none');
grid on
hold off
end
