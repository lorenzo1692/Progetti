function fig = plot_tf3d_global(gm, ttl, cpl)
%PLOT_TF3D_GLOBAL Results of the 3D beam + shell global model of a design point.
%
%   fig = PLOT_TF3D_GLOBAL(gm, ttl) plots, for coil 1, the output of
%   TF3D_GLOBAL_FROM_DESIGN:
%     (1) centreline coloured by the beam axial force and the deformed
%         shape (displacements magnified);
%     (2) axial force and in-plane bending moment along the coil;
%     (3) inner leg: vertical force of beam + vault shells against the
%         force T_bf of the 2D sizing;
%     (4) vault shells: hoop membrane force N11 and Tresca intensity.
%   fig = PLOT_TF3D_GLOBAL(gm, ttl, cpl) adds the history of the coupled
%   sizing (COUPLE_GLOBAL_SIZING) in a second figure.

if nargin < 2 || isempty(ttl), ttl = 'TF global model'; end
P = gm.P; M = size(P, 1);
U = gm.res.U(1:M, :);
Bf = gm.res.beam_f(1:M, :);
s = [0; cumsum(sqrt(sum(diff(P).^2, 2)))];
fig = figure('Name', [ttl ' - global model'], 'Color', 'w');

subplot(2, 2, 1); hold on; box on; axis equal
N = Bf(:, 1)/1e6;
scatter(P(:, 2), P(:, 3), 18, N, 'filled');
Uv = U(:, 1:3);
amp = 0.05*max(range(P(:, 2)), range(P(:, 3)))/max(max(abs(Uv(:))), eps);
plot([P(:, 2); P(1, 2)] + amp*[U(:, 2); U(1, 2)], [P(:, 3); P(1, 3)] + amp*[U(:, 3); U(1, 3)], 'k--');
cb = colorbar; ylabel(cb, 'axial force [MN]');
xlabel('R [m]'); ylabel('z [m]');
title(sprintf('coil 1: axial force; deformed x %.0f (max %.1f mm)', amp, 1e3*max(sqrt(sum(U(:, 1:3).^2, 2)))));

subplot(2, 2, 2); hold on; box on; grid on
plot(s, N, 'b-', 'LineWidth', 1.2);
plot(s, Bf(:, 6)/1e6, 'r-', 'LineWidth', 1.2);
xlabel('curvilinear abscissa from the inner mid plane [m]');
ylabel('MN, MN m');
legend({'beam axial force N (the OIS shells carry the rest)', 'in-plane moment M_z'}, 'Location', 'best');
title('beam forces along the coil');

subplot(2, 2, 3); hold on; box on; grid on
[z, o] = sort(gm.z_inner);
plot(gm.F_inner(o)/1e6, z, 'k-', 'LineWidth', 1.5);
plot(gm.N_beam(o)/1e6, z, 'b-');
plot(gm.N_vault(o)/1e6, z, 'm-');
yl = [min(z) max(z)];
plot(gm.T_bf/1e6*[1 1], yl, 'r--');
plot(gm.k_axial*gm.T_bf/1e6*[1 1], yl, 'r-');
xlabel('vertical force [MN]'); ylabel('z [m]');
legend({'inner leg total', 'beam', 'vault shells (one chord)', 'T_{bf} of the 2D sizing', ...
    sprintf('k T_{bf}, k = %.3f', gm.k_axial)}, 'Location', 'best');
title('straight inner leg: vertical force');

subplot(2, 2, 4); hold on; box on; grid on
plot(gm.vault_N11(o)/1e6, z, 'b-', 'LineWidth', 1.2);
plot(gm.vault_N22(o)/1e6, z, 'm-', 'LineWidth', 1.2);
xlabel('membrane force [MN/m]'); ylabel('z [m]');
legend({'N_{11} hoop', 'N_{22} vertical'}, 'Location', 'best');
title(sprintf('vault shells, t = %.0f mm, Tresca max %.0f MPa', 1e3*gm.section.t_vault, max(gm.vault_tresca)/1e6));

if nargin >= 3 && ~isempty(cpl) && numel(cpl.history) > 0
    h = cpl.history;
    f2 = figure('Name', [ttl ' - coupled sizing'], 'Color', 'w');
    it = 0:numel(h) - 1;
    subplot(1, 2, 1); hold on; box on; grid on
    plot(it, [h.k_sized], 'bo-'); plot(it, [h.k_global], 'rs-');
    xlabel('re-sizing'); ylabel('axial load factor k');
    legend({'k used by the 2D sizing', 'k from the global model'}, 'Location', 'best');
    title('coupled sizing: convergence');
    subplot(1, 2, 2); hold on; box on; grid on
    plot(it, 1e3*[h.radial_build], 'ko-'); plot(it, 1e3*[h.Nose], 'bs-');
    xlabel('re-sizing'); ylabel('mm');
    legend({'radial build', 'nose'}, 'Location', 'best');
    title('inner-leg build');
    fig = [fig f2];
end
end

function r = range(x)
r = max(x(:)) - min(x(:));
end
