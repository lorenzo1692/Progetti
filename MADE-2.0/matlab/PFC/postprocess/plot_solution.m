function plot_solution(DATA, sel_idx, tag, out_dir)
%PLOT_SOLUTION Plot and save the chosen PFC design point.
%
%   PLOT_SOLUTION(DATA, sel_idx, tag, out_dir) highlights design point
%   sel_idx over the full population (Iop vs WP_h0, colored by ST_MPa),
%   prints its full data row, and saves both a PNG of the highlighted
%   plot and an Excel file with that single row, named after tag (e.g.
%   the machine input file name and the PF coil number) and out_dir.
%
%   Mirrors MADE-2.0/matlab/CS/postprocess/plot_solution.m.
%
%   SCOPE NOTE (see manuale PFC): the legacy plot_PF_VNS.m also redrew the
%   full turn-by-turn cross-section (colored by |B| and by hoop stress)
%   for a single, hand-picked design point. That richer visualization was
%   not ported for v1 - deferred the same way CS's own plot_solution.m v1
%   deferred it (see manuale CS, "Fase 8").

if nargin < 4 || isempty(out_dir)
    out_dir = pwd;
end
if ~isfolder(out_dir)
    mkdir(out_dir);
end

row = DATA(sel_idx, :);
fprintf('\nSelected design point #%d (PF%g):\n', sel_idx, row.n_PF);
disp(row(:, intersect({'n_PF','Iop_kA','WP_h0','Btot','Ri','Re','WP_w', ...
    'n_layers','n_turns','SH_MPa','ST_MPa','L','E','Jeng'}, ...
    row.Properties.VariableNames, 'stable')));

fig = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
scatter(DATA.Iop_kA, DATA.WP_h0, 20, DATA.ST_MPa, 'filled');
hold on
scatter(row.Iop_kA, row.WP_h0, 250, 'kp', 'filled', 'MarkerEdgeColor', 'w', 'LineWidth', 1.5);
hold off
colormap(jet);
set(gca, 'FontSize', 16);
c = colorbar;
c.Label.String = 'Tresca stress S_T [MPa]';
xlabel('Iop [kA]', 'FontSize', 16);
ylabel('WP_h0 [m]', 'FontSize', 16);
title(sprintf('PF%g - selected design point #%d', row.n_PF, sel_idx), 'FontSize', 16);
grid on

fig_path = fullfile(out_dir, sprintf('%s_PF%g_design_%d.png', tag, row.n_PF, sel_idx));
print(fig, '-dpng', '-r300', fig_path);
fprintf('Plot saved to %s\n', fig_path);

xlsx_path = fullfile(out_dir, sprintf('%s_PF%g_design_%d.xlsx', tag, row.n_PF, sel_idx));
writetable(row, xlsx_path);
fprintf('Design point saved to %s\n', xlsx_path);
end
