function plot_solution(DATA, sel_idx, tag, out_dir)
%PLOT_SOLUTION Plot and save the chosen CS design point.
%
%   PLOT_SOLUTION(DATA, sel_idx, tag, out_dir) highlights design point
%   sel_idx over the full population (Iop vs Ri, colored by S_T_max),
%   prints its full data row, and saves both a PNG of the highlighted
%   plot and an Excel file with that single row, named after tag (e.g.
%   the machine input file name) and out_dir.
%
%   Mirrors MADE-2.0/matlab/postprocess/plot_solution.m (TF), with CS's
%   own natural coordinates (Ri instead of Rk_, S_T_max instead of
%   S_T_VT).
%
%   SCOPE NOTE (see manuale CS): the legacy Plot_CS_VNS.m also redrew the
%   full turn-by-turn module cross-section colored by local |B| (via a
%   second EMAG_FIELD_FORCES call). That richer visualization was not
%   ported for v1 - it can be added later as a separate function without
%   changing the solver, following TF's own plot_solution.m, which keeps
%   this same "highlight on the scan scatter" scope.

if nargin < 4 || isempty(out_dir)
    out_dir = pwd;
end
if ~isfolder(out_dir)
    mkdir(out_dir);
end

row = DATA(sel_idx, :);
fprintf('\nSelected design point #%d:\n', sel_idx);
disp(row(:, intersect({'Iop_kA','B_grades','Ri','Re','WP_w','WP_h', ...
    'n_layers','n_turns','S_hoop_max','S_T_max','L','E','Jeng'}, ...
    row.Properties.VariableNames, 'stable')));

fig = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
scatter(DATA.Iop_kA, DATA.Ri, 20, DATA.S_T_max, 'filled');
hold on
scatter(row.Iop_kA, row.Ri, 250, 'kp', 'filled', 'MarkerEdgeColor', 'w', 'LineWidth', 1.5);
hold off
colormap(jet);
set(gca, 'FontSize', 16);
c = colorbar;
c.Label.String = 'Tresca stress S_T [MPa]';
xlabel('Iop [kA]', 'FontSize', 16);
ylabel('Ri [m]', 'FontSize', 16);
title(sprintf('Selected design point #%d', sel_idx), 'FontSize', 16);
grid on

fig_path = fullfile(out_dir, sprintf('%s_design_%d.png', tag, sel_idx));
print(fig, '-dpng', '-r300', fig_path);
fprintf('Plot saved to %s\n', fig_path);

xlsx_path = fullfile(out_dir, sprintf('%s_design_%d.xlsx', tag, sel_idx));
writetable(row, xlsx_path);
fprintf('Design point saved to %s\n', xlsx_path);
end
