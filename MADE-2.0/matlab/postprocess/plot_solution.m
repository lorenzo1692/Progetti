function plot_solution(DATA, sel_idx, tag, out_dir)
%PLOT_SOLUTION Plot and save the chosen design point.
%
%   PLOT_SOLUTION(DATA, sel_idx, tag, out_dir) highlights design point
%   sel_idx over the full population (Iop vs Rk, colored by S_T_VT),
%   prints its full data row, and saves both a PNG of the highlighted
%   plot and an Excel file with that single row, named after tag (e.g.
%   the machine input file name) and out_dir.

if nargin < 4 || isempty(out_dir)
    out_dir = pwd;
end
if ~isfolder(out_dir)
    mkdir(out_dir);
end

row = DATA(sel_idx, :);
fprintf('\nSelected design point #%d:\n', sel_idx);
disp(row(:, intersect({'Iop','B_PHI_0','B_TF','Rk_','Rj_','Ri_','WP_h','WP_w', ...
    'n_layers','n_cond','S_T_VT','S_T_JT','L','E','Tau_discharge'}, ...
    row.Properties.VariableNames, 'stable')));

fig = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
scatter(DATA.Iop*1e-3, DATA.Rk_, 20, DATA.S_T_VT, 'filled');
hold on
scatter(row.Iop*1e-3, row.Rk_, 250, 'kp', 'filled', 'MarkerEdgeColor', 'w', 'LineWidth', 1.5);
hold off
colormap(jet);
set(gca, 'FontSize', 16);
c = colorbar;
c.Label.String = 'SINT VT [MPa]';
xlabel('Iop [kA]', 'FontSize', 16);
ylabel('Rk [m]', 'FontSize', 16);
title(sprintf('Selected design point #%d', sel_idx), 'FontSize', 16);
grid on

fig_path = fullfile(out_dir, sprintf('%s_design_%d.png', tag, sel_idx));
print(fig, '-dpng', '-r300', fig_path);

xlsx_path = fullfile(out_dir, sprintf('%s_design_%d.xlsx', tag, sel_idx));
writetable(row, xlsx_path);

fprintf('Saved selected design point #%d:\n  figure -> %s\n  data   -> %s\n', sel_idx, fig_path, xlsx_path);
end
