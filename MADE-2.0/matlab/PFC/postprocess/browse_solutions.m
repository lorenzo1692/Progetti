function sel_idx = browse_solutions(DATA)
%BROWSE_SOLUTIONS Show every feasible PFC design point (one coil, all WP_h) and let the user pick one.
%
%   sel_idx = BROWSE_SOLUTIONS(DATA) shows:
%     1) a summary uitable of all feasible design points found across the
%        WP_h sweep for one PF coil, and
%     2) an interactive scatter plot (Iop vs WP_h0, colored by ST_MPa)
%        where clicking a point preselects the closest design point.
%   It then asks in the terminal for the index of the design point to
%   plot and save (defaulting to the clicked point, if any), and returns
%   it as sel_idx.
%
%   Mirrors MADE-2.0/matlab/CS/postprocess/browse_solutions.m, with PFC's
%   own natural axes (WP_h0, the swept winding-pack height, instead of a
%   single geometric radius - see manuale PFC).

n = height(DATA);
if n == 0
    error('browse_solutions:no_solutions', 'DATA has no feasible design points to browse.');
end

fprintf('\n%d feasible design point(s) found.\n', n);

summary_cols = {'Index','Iop_kA','WP_h0','Btot','Ri','Re','WP_w','SH_MPa','ST_MPa','E'};
T_display = DATA(:, intersect(summary_cols, DATA.Properties.VariableNames, 'stable'));
T_display = addvars(T_display, (1:n)', 'Before', 1, 'NewVariableNames', 'Index');

f = uifigure('Name', 'WP PFC design points - summary', 'Position', [80 80 960 420]);
uitable(f, 'Data', T_display, 'Position', [10 10 940 400]);

fig = figure('Name', 'Click a point to preselect a design point (close the window to skip)', ...
    'units', 'normalized', 'outerposition', [0 0 0.6 0.6]); %#ok<NASGU>
scatter(DATA.Iop_kA, DATA.WP_h0, 40, DATA.ST_MPa, 'filled');
colormap(jet); colorbar;
xlabel('Iop [kA]'); ylabel('WP_h0 [m]');
title('Click a point to preselect a design (Tresca stress [MPa] as color)');
grid on

default_idx = 1;
try
    [xc, yc] = ginput(1);
    d2 = (DATA.Iop_kA - xc).^2 + (DATA.WP_h0 - yc).^2;
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
