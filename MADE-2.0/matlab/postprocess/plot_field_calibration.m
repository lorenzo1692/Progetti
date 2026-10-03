function fig = plot_field_calibration(cal, p, fig_title)
%PLOT_FIELD_CALIBRATION Start-of-scan peak-field calibration (WP_FIELD_CALIBRATION).
%
%   fig = PLOT_FIELD_CALIBRATION(cal, p, fig_title): left, the peak factor
%   k = peak on the conductor / Ampere field against the toroidal WP width
%   W, one line per operating current, with the constant corr_B_WP of the
%   smeared model for comparison (grid points without an admissible
%   reference layout, filled from the nearest one, are not drawn); right, the per-layer profile B/B_peak
%   against the fraction s of the turns from the layer inward, for the
%   narrowest and widest W at the middle current, against the linear
%   profile (B/B_peak = s) of the smeared model.

if nargin < 3, fig_title = ''; end
fig = figure('Name', 'Peak-field calibration', 'NumberTitle', 'off');
subplot(1, 2, 1); hold on; grid on
[W, o] = sort(cal.W);
cols = lines(numel(cal.Iop));
for j = 1:numel(cal.Iop)
    kj = cal.k(o, j); kj(cal.n_ref(o, j) == 0) = NaN;   % points filled from a neighbour: not drawn
    if all(isnan(kj)), continue, end
    plot(1e3*W, kj, '-o', 'Color', cols(j,:), 'LineWidth', 1.5, 'MarkerSize', 4, ...
        'DisplayName', sprintf('Iop = %.1f kA', cal.Iop(j)/1e3));
end
plot(1e3*[W(1) W(end)], p.corr_B_WP*[1 1], 'k--', 'LineWidth', 1.2, 'DisplayName', 'corr\_B\_WP (smeared)');
xlabel('Toroidal WP width W [mm]'); ylabel('k = B_{peak} / B_{Ampere}');
title('Peak factor'); legend('Location', 'best');
subplot(1, 2, 2); hold on; grid on
jm = ceil(numel(cal.Iop)/2);
plot(cal.s, squeeze(cal.f(o(1), jm, :)), 'LineWidth', 1.5, 'DisplayName', sprintf('W = %.0f mm', 1e3*W(1)));
plot(cal.s, squeeze(cal.f(o(end), jm, :)), 'LineWidth', 1.5, 'DisplayName', sprintf('W = %.0f mm', 1e3*W(end)));
plot([0 1], [0 1], 'k--', 'LineWidth', 1.2, 'DisplayName', 'linear (smeared)');
xlabel('s = turns from the layer inward / total'); ylabel('B_{layer} / B_{peak}');
title(sprintf('Layer profile, Iop = %.1f kA', cal.Iop(jm)/1e3)); legend('Location', 'southeast');
if ~isempty(fig_title)
    annotation(fig, 'textbox', [0 0.94 1 0.06], 'String', [fig_title ' - peak-field calibration'], ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
end
end
