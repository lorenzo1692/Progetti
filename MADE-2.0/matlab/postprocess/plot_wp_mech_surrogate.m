function plot_wp_mech_surrogate(out, p, fig_title)
%PLOT_WP_MECH_SURROGATE Stress map and figure of merit of WP_MECH_SURROGATE.
%
%   PLOT_WP_MECH_SURROGATE(out, p) plots, for one solved design point:
%     1. the Tresca stress intensity over the whole inner-leg section
%        (element centres, MPa) with the case stress-classification lines;
%     2. the jacket stresses per layer (maximum over the turns): primary
%        membrane Pm and membrane+bending Pm+Pb (Lorentz + axial load case)
%        against Sm and 1.5 Sm, Pm+Pb under all loads (P+Q, with the
%        cool-down) and the peak in the cable fillet (for information);
%   and prints the validity checks and the figure of merit in the console.
%   A result whose checks fail is labelled INVALID.
%
%   out - result of WP_MECH_SURROGATE; p - machine parameters.

if nargin < 3 || isempty(fig_title), fig_title = 'WP mechanical surrogate'; end
Smj = out.fom.Sm_jacket; Smc = out.fom.Sm_case;
m = out.mesh; xy = m.xy;

figure('units', 'normalized', 'outerposition', [0.05 0.08 0.9 0.84]);

% --- 1. Tresca map: whole section and WP zoom ---------------------------
clim_ = [0, max(1.5*max(Smj, Smc), out.fom.jacket_peak)/1e6];
ax1 = subplot(1, 3, 1);
draw_map(ax1, out, clim_);
for q = 1:numel(out.case_scl)
    s = out.case_scl(q);
    plot(ax1, [s.P0(1) s.P1(1)], [s.P0(2) s.P1(2)], 'k-', 'LineWidth', 1.2);
    text(ax1, s.P1(1), s.P1(2), sprintf(' %.0f/%.0f', s.Pm/1e6, s.PmPb/1e6), 'FontSize', 7);
end
title(ax1, {fig_title, 'steel Tresca [MPa]; case SCLs: Pm/Pm+Pb'}, 'Interpreter', 'none');

ax2 = subplot(1, 3, 2);
draw_map(ax2, out, clim_);
[~, it] = max([out.turn.peak]);
plot(ax2, out.turn(it).peak_xy(1), out.turn(it).peak_xy(2), 'kp', 'MarkerFaceColor', 'w', 'MarkerSize', 12);
cav = out.geo.cav;
xlim(ax2, [min(cav(:,1)) max(cav(:,1))] + [-0.01 0.01]);
ylim(ax2, [min(cav(:,2)) max(cav(:,2))] + [-0.01 0.01]);
cb = colorbar(ax2); ylabel(cb, 'Tresca stress intensity [MPa]');
title(ax2, sprintf('WP: jacket peak %.0f MPa (star)', out.fom.jacket_peak/1e6));

% --- 2. Jacket per layer ------------------------------------------------
f = out.fom;
Smj = f.Sm_jacket; Smc = f.Sm_case;
if f.classified, PL = out.primary.layer; else, PL = out.layer; end
ax3 = subplot(1, 3, 3);
nl = numel(PL.Pm);
hb = bar(ax3, 1:nl, [PL.Pm; PL.PmPb]'/1e6, 'grouped');
hold(ax3, 'on');
hq = plot(ax3, 1:nl, out.layer.PmPb/1e6, 'ks--', 'MarkerFaceColor', 'w');
hp = plot(ax3, 1:nl, out.layer.peak/1e6, 'kd-', 'MarkerFaceColor', 'k');
hs = plot(ax3, [0.5 nl+0.5], [Smj Smj]/1e6, 'r--', 'LineWidth', 1.2);
h15 = plot(ax3, [0.5 nl+0.5], 1.5*[Smj Smj]/1e6, 'r:', 'LineWidth', 1.2);
hold(ax3, 'off');
xlabel(ax3, 'Layer (1 = plasma side)'); ylabel(ax3, 'Jacket stress [MPa]');
legend(ax3, [hb(1) hb(2) hq hp hs h15], {'P_m primary', 'P_m+P_b primary', ...
    'P_m+P_b all loads (P+Q)', 'peak in cable fillet', 'S_m', '1.5 S_m'}, 'Location', 'northwest');
title(ax3, {sprintf('Jacket P_m %.0f, P_m+P_b %.0f, P+Q %.0f, peak %.0f MPa', f.jacket_Pm/1e6, ...
    f.jacket_PmPb/1e6, f.jacket_PQ/1e6, f.jacket_peak/1e6), ...
    sprintf('Case P_m %.0f, P_m+P_b %.0f MPa - %s', f.case_Pm/1e6, f.case_PmPb/1e6, f.status)}, ...
    'Interpreter', 'tex');
grid(ax3, 'on');

% --- console table --------------------------------------------------
fprintf('\nMechanical surrogate - checks: %s\n', out.checks.summary);
c = out.checks;
fprintf(['  area %.1e | min detJ %.2e | contact converged %d | residual %.1e | ' ...
    'force balance %.1e | axial %.1e | SCL coverage %.2f\n'], c.area_rel_err, c.min_detJ, ...
    c.contact_converged, c.solve_residual, c.force_balance, c.axial_rel_err, c.scl_min_coverage);
if f.classified
    fprintf('Figure of merit (linearized Tresca; P = Lorentz + axial, Q = cool-down):\n');
else
    fprintf('Figure of merit (UNCLASSIFIED: all loads used as primary, conservative):\n');
end
for k = 1:numel(f.util)
    fprintf('  %-24s %.3f\n', f.util_names{k}, f.util(k));
end
fprintf('  jacket Pm    %5.0f MPa at layer %d col %d (Sm = %.0f MPa)\n', f.jacket_Pm/1e6, f.jacket_Pm_at, Smj/1e6);
fprintf('  jacket Pm+Pb %5.0f MPa at layer %d col %d (1.5 Sm = %.0f MPa)\n', f.jacket_PmPb/1e6, f.jacket_PmPb_at, 1.5*Smj/1e6);
fprintf('  jacket P+Q   %5.0f MPa at layer %d col %d (3 Sm = %.0f MPa)\n', f.jacket_PQ/1e6, f.jacket_PQ_at, 3*Smj/1e6);
fprintf('  jacket peak  %5.0f MPa at layer %d col %d (local, for information)\n', f.jacket_peak/1e6, f.jacket_peak_at);
fprintf('  case Pm %.0f MPa (%s), Pm+Pb %.0f MPa (%s), P+Q %.0f MPa (%s); Sm = %.0f MPa\n', ...
    f.case_Pm/1e6, f.case_Pm_at, f.case_PmPb/1e6, f.case_PmPb_at, f.case_PQ/1e6, f.case_PQ_at, Smc/1e6);
fprintf('  -> %s\n', f.status);
end

function draw_map(ax, out, clim_)
% steel (jacket, case) coloured by Tresca, other materials in grey
m = out.mesh; xy = m.xy;
hold(ax, 'on');
qs = m.q8_mat == 2; ts = m.t6_mat == 4;
grey = [0.82 0.82 0.82];
patch(ax, 'Faces', m.q8(~qs, [1 5 2 6 3 7 4 8]), 'Vertices', xy, 'FaceColor', grey, 'EdgeColor', 'none');
patch(ax, 'Faces', m.t6(~ts, [1 4 2 5 3 6]), 'Vertices', xy, 'FaceColor', grey, 'EdgeColor', 'none');
patch(ax, 'Faces', m.q8(qs, [1 5 2 6 3 7 4 8]), 'Vertices', xy, ...
    'FaceVertexCData', out.elem_SINT.q8(qs)/1e6, 'FaceColor', 'flat', 'EdgeColor', 'none');
patch(ax, 'Faces', m.t6(ts, [1 4 2 5 3 6]), 'Vertices', xy, ...
    'FaceVertexCData', out.elem_SINT.t6(ts)/1e6, 'FaceColor', 'flat', 'EdgeColor', 'none');
axis(ax, 'equal'); box(ax, 'on');
colormap(ax, jet(256));
set(ax, 'CLim', clim_);
xlabel(ax, 'Toroidal x [m]'); ylabel(ax, 'Radial y [m]');
end
