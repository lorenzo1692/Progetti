function grades = plot_hotspot_transient(row, p, fig_title)
%PLOT_HOTSPOT_TRANSIENT Hot-spot temperature T(t) after a quench, one curve per grade.
%
%   grades = PLOT_HOTSPOT_TRANSIENT(row, p) re-runs, for every conductor
%   grade of the design point, the same adiabatic power-balance ODE that
%   CICC used to size its copper (HEAT_BALANCE_CICC_ODE, constants from
%   CICC_PARAMS) and plots the full temperature history T(t): flat-top
%   current for Tau_delay, then exponential discharge with Tau_discharge,
%   with T saturating at the hot-spot value. One line per grade, the
%   type-specific allowable (p.THS_max_LTS / p.THS_max_HTS) as dashed
%   lines, and per grade the same summary the original heat-balance plot
%   printed (THS, Cu areas, SC area, Nsc/NCu, Tau_del/Tau_dis, Iop/Bop).
%
%   A grade is a run of consecutive layers with the same cable
%   (type_cable, N_Sc, N_Cu). Its operating field is row.B_grade (the
%   field CICC sized it at, stored by SCAN_WP_DESIGNS); for rows without
%   B_grade, the smeared B_layers value at the grade's first layer is used.
%
%   row - DATA row or manual struct. Used: Iop, n_layers, n_turns,
%         type_cable, N_Sc, N_Cu, Tau_discharge, B_TF (+ B_grade if present).
%   p   - machine parameter struct (THS_max_LTS, THS_max_HTS).
%
%   grades - struct array: layers, type, mat, B, N_Sc, N_Cu, THS, t, T.

if nargin < 3 || isempty(fig_title)
    fig_title = sprintf('Hot-spot transient - Iop=%.0f A, Tau_{dis}=%.1f s', row.Iop, row.Tau_discharge);
end

cp = cicc_params();
n_layers = row.n_layers;
n_turns = row.n_turns(1:n_layers);
type_cable = row.type_cable(1:n_layers);
N_Sc = row.N_Sc(1:n_layers);
N_Cu = row.N_Cu(1:n_layers);

% Field each layer's cable was sized at
if isfield_or_var(row, 'B_grade') && any(row.B_grade(1:n_layers) > 0)
    B_lay = row.B_grade(1:n_layers);
else
    n_spire_ = zeros(1, n_layers);
    n_spire_(1) = sum(n_turns);
    for k = 2:n_layers
        n_spire_(k) = n_spire_(k-1) - n_turns(k);
    end
    B_lay = row.B_TF .* (n_spire_/n_spire_(1));
end

% Split into grades (runs of identical cable)
starts = 1;
for k = 2:n_layers
    if ~strcmp(type_cable{k}, type_cable{k-1}) || N_Sc(k) ~= N_Sc(k-1) || N_Cu(k) ~= N_Cu(k-1)
        starts(end+1) = k; %#ok<AGROW>
    end
end
ends = [starts(2:end)-1, n_layers];

n_g = numel(starts);
grades = struct('layers', {}, 'type', {}, 'mat', {}, 'B', {}, 'N_Sc', {}, 'N_Cu', {}, ...
    'THS', {}, 't', {}, 'T', {});
for g = 1:n_g
    k = starts(g);
    B = B_lay(k);
    if strcmp(type_cable{k}, 'HTS')
        mat = 2; sc_name = 'REBCO'; lim = p.THS_max_HTS;
    elseif B < cp.B_NbTi_max
        mat = 1; sc_name = 'NbTi';  lim = p.THS_max_LTS;
    else
        mat = 0; sc_name = 'Nb_3Sn'; lim = p.THS_max_LTS;
    end
    [THS, t, T] = heat_balance_cicc_ode(N_Sc(k), N_Cu(k), cp.d_fili, cp.CunonCu, row.Iop, B, ...
        row.Tau_discharge, mat, cp.d_cc, cp.VF, cp.cos_theta, cp.S_tapes, cp.Tau_delay);
    grades(g).layers = starts(g):ends(g);
    grades(g).type = type_cable{k};
    grades(g).sc_name = sc_name;
    grades(g).mat = mat;
    grades(g).B = B;
    grades(g).N_Sc = N_Sc(k);
    grades(g).N_Cu = N_Cu(k);
    grades(g).THS = THS;
    grades(g).limit = lim;
    grades(g).t = t;
    grades(g).T = T;
end

%% Plot
figure('units', 'normalized', 'outerposition', [0.1 0.1 0.75 0.75]);
ax = axes('Position', [0.07 0.10 0.55 0.80]);
hold(ax, 'on');
cols = lines(max(n_g, 1));
h = gobjects(1, n_g);
labels = cell(1, n_g);
t_end = 0;
for g = 1:n_g
    G = grades(g);
    labels{g} = sprintf('Grade %d: L%d-%d %s %s, B=%.2f T, THS=%.0f K', g, ...
        G.layers(1), G.layers(end), G.type, strrep(G.sc_name, '_', ''), G.B, G.THS);
    h(g) = plot(ax, G.t, G.T, '-', 'Color', cols(g,:), 'LineWidth', 1.8);
    [~, im] = max(G.T);
    plot(ax, G.t(im), G.T(im), 'o', 'Color', cols(g,:), 'MarkerFaceColor', cols(g,:), ...
        'HandleVisibility', 'off');
    t_end = max(t_end, G.t(end));
end
% Allowables of the cable types present, and the discharge start
lims = unique([grades.limit]);
for L = lims
    types = unique({grades([grades.limit] == L).type});
    plot(ax, [0 t_end], [L L], 'r--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
    text(ax, 0.02*t_end, L, sprintf(' %s limit %.0f K', strjoin(types, '/'), L), ...
        'Color', 'r', 'VerticalAlignment', 'bottom', 'FontSize', 9);
end
plot(ax, [cp.Tau_delay cp.Tau_delay], [0 max([lims, grades.THS])*1.05], 'k:', 'HandleVisibility', 'off');
text(ax, cp.Tau_delay, 5, ' discharge start', 'FontSize', 8, 'VerticalAlignment', 'bottom');
hold(ax, 'off');
xlim(ax, [0 t_end]);
ylim(ax, [0 max([lims, grades.THS])*1.1]);
xlabel(ax, 't [s]'); ylabel(ax, 'T_{hot spot} [K]');
title(ax, fig_title);
grid(ax, 'on');
legend(h, labels, 'Location', 'southeast', 'FontSize', 8);

% Per-grade summary, as in the original heat-balance plot
ax2 = axes('Position', [0.65 0.05 0.34 0.90], 'Visible', 'off');
y = 1;
dy = 0.95/max(n_g, 1);
for g = 1:n_g
    G = grades(g);
    if G.mat == 2
        A_sc = G.N_Sc*cp.S_tapes;
        A_cu_ns = 0;
    else
        A_sc = G.N_Sc*pi*cp.d_fili^2/(4*(1+cp.CunonCu));
        A_cu_ns = A_sc;
    end
    A_cu_s = G.N_Cu*pi*cp.d_fili^2/4;
    keyt = {sprintf('\\bfGrade %d - layers %d-%d (%s)\\rm', g, G.layers(1), G.layers(end), G.type), ...
        sprintf('Hot spot Temperature = %.0f [K]  (limit %.0f K)', G.THS, G.limit), ...
        sprintf('Cu Area non seg = %.1f [mm^2]', A_cu_ns*1e6), ...
        sprintf('Cu Area seg = %.1f [mm^2]', A_cu_s*1e6), ...
        sprintf('%s Area = %.1f [mm^2]', G.sc_name, A_sc*1e6), ...
        sprintf('N_{sc} = %d   N_{Cu} = %d', G.N_Sc, G.N_Cu), ...
        sprintf('\\tau_{del} = %.1f [s]   \\tau_{dis} = %.2f [s]', cp.Tau_delay, row.Tau_discharge), ...
        sprintf('I_{op} = %.2f [kA]   B_{op} = %.2f [T]', row.Iop/1000, G.B)};
    text(ax2, 0, y, keyt, 'Units', 'normalized', 'VerticalAlignment', 'top', ...
        'FontSize', 9, 'Color', cols(g,:)*0.8);
    y = y - dy;
end

fprintf('Hot-spot transient:\n');
for g = 1:n_g
    fprintf('  grade %d (L%d-%d, %s %s, B=%.2f T): THS = %.0f K (limit %.0f K), N_Sc=%d, N_Cu=%d\n', ...
        g, grades(g).layers(1), grades(g).layers(end), grades(g).type, grades(g).sc_name, ...
        grades(g).B, grades(g).THS, grades(g).limit, grades(g).N_Sc, grades(g).N_Cu);
end
end

function tf = isfield_or_var(row, name)
% true if a struct field or table variable called NAME exists
if istable(row)
    tf = any(strcmp(row.Properties.VariableNames, name));
else
    tf = isfield(row, name);
end
end
