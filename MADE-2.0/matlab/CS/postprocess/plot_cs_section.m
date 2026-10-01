function plot_cs_section(DATA, sel_idx, tag, p, g, out_dir)
%PLOT_CS_SECTION Plot one module and the full CS stack, colored by local |B|.
%
%   PLOT_CS_SECTION(DATA, sel_idx, tag, p, g, out_dir) redraws the chosen
%   design point turn by turn: every cable is drawn as turn insulation,
%   jacket and conductor core, with the core colored by the local field
%   |B| = sqrt(Br^2+Bz^2) from EMAG_FIELD_FORCES (full_stack=1). Two PNGs
%   are saved: "<tag>_design_<idx>_module.png" (one module, with
%   dimensions and key numbers) and "<tag>_design_<idx>_stack.png" (all
%   n_moduli modules). Both share the same color scale.
%
%   Ported from the plotting part of the legacy plot_CS_VNS.m. DATA does
%   not store the per-grade radii, so they are rebuilt with the same
%   recurrence used by SCAN_WP_DESIGNS (Ri_grades(v) = Re_grades(v) -
%   Cond_w(v)*n_layers(v), Re_grades(v+1) = Ri_grades(v) - tins), which is
%   why p and g are needed. Only the rectangular cable (shape_cable=201)
%   is drawn, as in the scan.

if nargin < 6 || isempty(out_dir)
    out_dir = pwd;
end
if ~isfolder(out_dir)
    mkdir(out_dir);
end

row = DATA(sel_idx, :);
n_grades = p.n_grades;
n_moduli = p.n_moduli;
n_layers = row.n_layers(1, :);
n_turns  = row.n_turns(1, :);
Cond_w = row.Cond_w(1, :);
Cond_h = row.Cond_h(1, :);
JT = row.JT(1, :);
Iop = row.Iop_kA * 1e3;
tins = g.tins;
grins_h = g.grins_h;
WP_h = g.WP_h;
Ri = row.Ri;
Re = row.Re;

Re_grades = zeros(1, n_grades+1);
Ri_grades = zeros(1, n_grades);
Re_grades(1) = g.Re_WP_outer;
for v = 1:n_grades
    Ri_grades(v) = Re_grades(v) - Cond_w(v)*n_layers(v);
    Re_grades(v+1) = Ri_grades(v) - tins;
end

[~, ~, ~, ~, ~, BR, BZ, ~, ~, ~, ~] = emag_field_forces( ...
    Cond_h, Cond_w, Ri_grades, n_turns, n_layers, n_grades, n_moduli, Iop, g.spacer, 1);
B_turn = sqrt(BR.^2 + BZ.^2);
B_lo = min(B_turn);
B_hi = max(B_turn);
n_col = 256;
map = jet(n_col);
if B_hi > B_lo
    idx_of = @(b) min(n_col, max(1, round(1 + (b - B_lo)/(B_hi - B_lo)*(n_col-1))));
else
    idx_of = @(b) 1;
end

n_per_module = numel(B_turn) / n_moduli;

%% Figure 1: one module
mod_plot = min(3, n_moduli);
fig1 = figure('units', 'normalized', 'outerposition', [0.25 0.1 0.5 0.9]);
draw_module(0, mod_plot);
axis equal
format_axes(B_lo, B_hi, map);
title(sprintf('CS design #%d - module %d - %s', sel_idx, mod_plot, tag), 'Interpreter', 'none', 'FontSize', 14);

plot([0 Ri], [0 0], '--.k', 'linewidth', 1);
plot([0 Re], [-0.05 -0.05], '--.k', 'linewidth', 1);
plot([Ri Re], [WP_h+0.01 WP_h+0.01], '--.b', 'linewidth', 1);
plot([Re+0.01 Re+0.01], [0 WP_h], '--.b', 'linewidth', 1);
text(Ri, -0.03, sprintf('Ri=%.3g m', Ri), 'HorizontalAlignment', 'center', 'FontSize', 12);
text(Re, -0.08, sprintf('Re=%.3g m', Re), 'HorizontalAlignment', 'center', 'FontSize', 12);
text(Re+0.02, WP_h/2, sprintf('h=%.3g m', WP_h), 'FontSize', 12);
text((Ri+Re)/2, WP_h+0.03, sprintf('w=%.3g m', row.WP_w), 'HorizontalAlignment', 'center', 'FontSize', 12);

keyt = { ...
    sprintf('\\sigma_{hoop} = %.1f MPa', row.S_hoop_max), ...
    sprintf('\\sigma_{T} = %.1f MPa', row.S_T_max), ...
    sprintf('\\Phi = %.1f Wb', row.Phi_TOT), ...
    sprintf('Layers: %s', mat2str(n_layers)), ...
    sprintf('Turns: %s', mat2str(n_turns)), ...
    sprintf('Iop = %.2f kA', row.Iop_kA), ...
    sprintf('Jeng = %.1f A/mm^2', row.Jeng), ...
    sprintf('L = %.4g H', row.L), ...
    sprintf('E = %.1f MJ', row.E), ...
    sprintf('dx = %s mm', mat2str(round(Cond_w*1e3, 1))), ...
    sprintf('dy = %s mm', mat2str(round(Cond_h*1e3, 1))), ...
    sprintf('jt = %s mm', mat2str(round(JT*1e3, 1))), ...
    sprintf('Fz = %.1f MN', row.Fz_MN), ...
    sprintf('V = %.1f kV', p.V_MAX/1e3)};
text(0.01, 0.6*WP_h, keyt, 'HorizontalAlignment', 'left', 'FontSize', 11, 'Interpreter', 'tex');
xlim([0, Re+0.25]);
ylim([-0.25, WP_h+0.25]);

f1 = fullfile(out_dir, sprintf('%s_design_%d_module.png', tag, sel_idx));
print(fig1, '-dpng', '-r300', f1);
fprintf('Module section saved to %s\n', f1);

%% Figure 2: full stack
fig2 = figure('units', 'normalized', 'outerposition', [0.25 0.1 0.5 0.9]);
for m = 1:n_moduli
    draw_module((WP_h + g.spacer)*(m-1), m);
end
axis equal
format_axes(B_lo, B_hi, map);
title(sprintf('CS design #%d - full stack - %s', sel_idx, tag), 'Interpreter', 'none', 'FontSize', 14);
f2 = fullfile(out_dir, sprintf('%s_design_%d_stack.png', tag, sel_idx));
print(fig2, '-dpng', '-r300', f2);
fprintf('Full CS section saved to %s\n', f2);

    function draw_module(dz, m)
        plot([Ri Re Re Ri Ri], [dz dz dz+WP_h dz+WP_h dz], 'k', 'linewidth', 1); hold on
        n = (m-1)*n_per_module;
        for v = 1:n_grades
            for t = 1:n_turns(v)
                for l = 1:n_layers(v)
                    n = n + 1;
                    x0 = Ri_grades(v) + Cond_w(v)*(l-1);
                    y0 = Cond_h(v)*(t-1) + grins_h + dz;
                    rectangle('Position', [x0, y0, Cond_w(v), Cond_h(v)], 'FaceColor', 'w');
                    rectangle('Position', [x0+tins, y0+tins, Cond_w(v)-2*tins, Cond_h(v)-2*tins], ...
                        'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'w');
                    s = JT(v) + tins;
                    rectangle('Position', [x0+s, y0+s, Cond_w(v)-2*s, Cond_h(v)-2*s], ...
                        'FaceColor', map(idx_of(B_turn(n)), :));
                end
            end
        end
    end
end

function format_axes(B_lo, B_hi, map)
set(gca, 'FontSize', 14);
colormap(map);
c = colorbar;
c.Label.String = '|B| [T]';
c.Label.FontSize = 14;
if B_hi > B_lo
    clim([B_lo B_hi]);
end
xlabel('r [m]'); ylabel('z [m]');
end
