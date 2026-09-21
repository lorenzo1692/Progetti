function sel_idx = browse_solutions(DATA)
%BROWSE_SOLUTIONS Show every feasible design point and let the user pick one.
%
%   sel_idx = BROWSE_SOLUTIONS(DATA) shows:
%     1) a summary uitable of all feasible design points found by
%        SCAN_WP_DESIGNS, and
%     2) an interactive scatter plot (Iop vs Rk, colored by S_T_VT) where
%        clicking a point preselects the closest design point.
%   It then asks in the terminal for the index of the design point to
%   plot and save (defaulting to the clicked point, if any), and returns
%   it as sel_idx.

n = height(DATA);
if n == 0
    error('browse_solutions:no_solutions', 'DATA has no feasible design points to browse.');
end

fprintf('\n%d feasible design point(s) found.\n', n);

summary_cols = {'Index','Iop','B_PHI_0','Rk_','WP_h','WP_w','n_layers','n_cond','S_T_VT','S_T_JT'};
T_display = DATA(:, intersect(summary_cols, DATA.Properties.VariableNames, 'stable'));
T_display = addvars(T_display, (1:n)', 'Before', 1, 'NewVariableNames', 'Index');

f = uifigure('Name', 'WP TF design points - summary', 'Position', [80 80 920 420]);
uitable(f, 'Data', T_display, 'Position', [10 10 900 400]);

fig = figure('Name', 'Click a point to preselect a design point (close the window to skip)', ...
    'units', 'normalized', 'outerposition', [0 0 0.6 0.6]);
scatter(DATA.Iop*1e-3, DATA.Rk_, 40, DATA.S_T_VT, 'filled');
colormap(jet); colorbar;
xlabel('Iop [kA]'); ylabel('Rk [m]');
title('Click a point to preselect a design (Tresca stress on the Jacket [MPa] as color)');
grid on

default_idx = 1;
try
    [xc, yc] = ginput(1);
    d2 = (DATA.Iop*1e-3 - xc).^2 + (DATA.Rk_ - yc).^2;
    [~, default_idx] = min(d2);
    fprintf('Closest design point to your click: #%d\n', default_idx);
catch
    fprintf('No click registered; defaulting to design point #%d.\n', default_idx);
end

answer = input(sprintf('Enter the index of the design point to plot and save [%d]: ', default_idx), 's');
if isempty(answer)
    sel_idx = default_idx;
else
    sel_idx = round(str2double(answer));
end

if isnan(sel_idx) || sel_idx < 1 || sel_idx > n
    error('browse_solutions:invalid_index', 'Invalid index: must be an integer between 1 and %d.', n);
end

fprintf('Selected design point #%d.\n', sel_idx);
end
