function plot_tf3d(c3, p, fig_title)
%PLOT_TF3D Plot the 3D coil shape and field checks from TF3D_FROM_DESIGN.
%
%   PLOT_TF3D(c3, p, fig_title) draws four panels: the bending-free shape
%   with its three-arc fit, a 3D view of one coil with the |B| color scale
%   from the single-filament model, the toroidal field ripple at the
%   outboard midplane, and a text summary of the checks against MADE's own
%   2D quantities (Ampere's law, T_bf, the field peak at the inner-leg
%   cross-section).

if nargin < 3 || isempty(fig_title), fig_title = '3D coil shape and field'; end

figure('units', 'normalized', 'outerposition', [0.05 0.05 0.9 0.85]);

%% Panel 1: shape + arc fit
subplot(2,2,1)
plot(c3.shape.r, c3.shape.z, 'k-', 'LineWidth', 1.2); hold on
plot(c3.shape.r, -c3.shape.z, 'k-', 'LineWidth', 1.2);
th1 = linspace(c3.arcs.ang1(1), c3.arcs.ang1(2), 40);
th2 = linspace(c3.arcs.ang2(1), c3.arcs.ang2(2), 40);
th3 = linspace(c3.arcs.ang3(1), c3.arcs.ang3(2), 40);
plot(c3.arcs.C1(1)+c3.arcs.R1*cos(th1), c3.arcs.C1(2)+c3.arcs.R1*sin(th1), 'r--', 'LineWidth', 1.5);
plot(c3.arcs.C2(1)+c3.arcs.R2*cos(th2), c3.arcs.C2(2)+c3.arcs.R2*sin(th2), 'b--', 'LineWidth', 1.5);
plot(c3.arcs.C3(1)+c3.arcs.R3*cos(th3), c3.arcs.C3(2)+c3.arcs.R3*sin(th3), 'g--', 'LineWidth', 1.5);
plot(c3.arcs.C1(1)+c3.arcs.R1*cos(th1), -(c3.arcs.C1(2)+c3.arcs.R1*sin(th1)), 'r--', 'LineWidth', 1);
plot(c3.arcs.C2(1)+c3.arcs.R2*cos(th2), -(c3.arcs.C2(2)+c3.arcs.R2*sin(th2)), 'b--', 'LineWidth', 1);
plot(c3.arcs.C3(1)+c3.arcs.R3*cos(th3), -(c3.arcs.C3(2)+c3.arcs.R3*sin(th3)), 'g--', 'LineWidth', 1);
plot(c3.geo.r1, c3.shape.z(1), 'ks', 'MarkerFaceColor', 'k');
plot(c3.geo.r2, 0, 'ks', 'MarkerFaceColor', 'k');
axis equal; grid on
xlabel('r [m]'); ylabel('z [m]');
title(sprintf('Shape mode %d - arc fit RMS %.1f mm, max %.1f mm', c3.shape.mode, 1e3*c3.arcs.rms, 1e3*c3.arcs.max_dev));
legend('bending-free', '', 'arc1', 'arc2', 'arc3', 'Location', 'best');

%% Panel 2: 3D view
subplot(2,2,2)
n_TF = c3.geo.n_TF;
hold on
for c = 1:n_TF
    a = (c-1)*2*pi/n_TF;
    ca = cos(a); sa = sin(a);
    rr = c3.path_rz(:,1); zz = c3.path_rz(:,2);
    x = [rr; flipud(rr)]; z = [zz; -flipud(zz)];
    X = x*ca; Y = x*sa; Z = z;
    plot3(X, Y, Z, 'b-', 'LineWidth', 1);
end
axis equal; grid on; view(35, 20)
xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
title(sprintf('%d TF coils - current centroid line', n_TF));

%% Panel 3: ripple bar
subplot(2,2,3)
bar([c3.checks.ripple, c3.checks.ripple_target]);
set(gca, 'XTickLabel', {'3D model', 'target (p.ripple)'});
ylabel('ripple at outboard midplane');
grid on
title('Field ripple check');

%% Panel 4: text summary
subplot(2,2,4); axis off
txt = {
    fig_title, '', ...
    sprintf('Ampere check at R0: B ratio = %.4f (target 1.0000)', c3.checks.ampere_ratio), ...
    sprintf('Ripple: %.4f (target %.4f)', c3.checks.ripple, c3.checks.ripple_target), ...
    sprintf('T axial / T_bf (2D shell model): %.4f', c3.checks.T_ratio), ...
    sprintf('Peak field, inner-leg x-section, 3D/2D: %.4f', c3.checks.peak_ratio_mid), ...
    sprintf('  2D peak (validated): %.3f T', c3.checks.peak2d_mid), ...
    sprintf('  3D-corrected estimate: %.3f T', c3.checks.peak3d_mid), ...
    sprintf('Arc fit: RMS %.1f mm, max %.1f mm', 1e3*c3.arcs.rms, 1e3*c3.arcs.max_dev)
    };
text(0.02, 0.95, txt, 'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 11, 'Interpreter', 'none');

fprintf('\n%s\n', c3.checks.summary);
end
