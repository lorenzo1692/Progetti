function fem_plot_results(res, icase, out_png)
%FEM_PLOT_RESULTS Field and stress maps of one FEM load case, with the design-model comparison.
%
%   FEM_PLOT_RESULTS(res, icase, out_png) draws, for load case icase of
%   res (from FEM_PFC_VERIFY; default 1 = design case), the winding-pack
%   cross-section colored by the vertical field Bz, the jacket hoop
%   stress, the jacket vertical stress and the Tresca stress, plus the
%   hoop stress per layer against the analytic ring model. Saves a PNG if
%   out_png is given.

if nargin < 2 || isempty(icase), icase = 1; end
c = res.cases(icase);
mesh = res.mesh; wp = res.wp; ne = size(mesh.elems, 1);
Bz_el = mean(c.BZ, 2);
S_T = zeros(ne, 1);
for e = 1:ne
    sg = [c.s_th(e), c.s_z(e), c.s_r(e)];
    S_T(e) = max(sg) - min(sg);
end

fig = figure('Name', sprintf('FEM PF%d - %s', wp.n_PF, c.name), 'units', 'normalized', 'outerposition', [0.05 0.05 0.9 0.85]);
maps = {Bz_el, c.s_th*1e-6, c.s_z*1e-6, S_T*1e-6};
titles = {'B_z [T]', 'Jacket hoop stress [MPa]', 'Jacket vertical stress [MPa]', 'Tresca [MPa]'};
for k = 1:4
    subplot(2, 3, k)
    patch('Faces', mesh.elems, 'Vertices', mesh.nodes, 'FaceVertexCData', maps{k}, 'FaceColor', 'flat', 'EdgeColor', 'none');
    axis equal tight; colormap(jet); colorbar;
    title(titles{k}); xlabel('r [m]'); ylabel('z (local) [m]');
end

subplot(2, 3, 5)
plot(1:wp.n_l, c.hoop_layer_MPa, 'o-', 'LineWidth', 1.5); hold on
if icase == 1
    plot(1:wp.n_l, res.analytic.hoop_layer_MPa, 's--', 'LineWidth', 1.5);
    legend('FEM', 'analytic ring', 'Location', 'best');
end
grid on; xlabel('layer (1 = inner)'); ylabel('hoop stress [MPa]'); title('Hoop stress per layer');

subplot(2, 3, 6)
names = {res.cases.name};
vals = [res.cases.hoop_max_MPa; res.cases.S_T_tresca_MPa];
bar(vals');
set(gca, 'XTick', 1:numel(names), 'XTickLabel', names); xtickangle(45);
legend('hoop', 'Tresca', 'Location', 'best'); grid on; ylabel('MPa'); title('All load cases');

if nargin >= 3 && ~isempty(out_png)
    print(fig, '-dpng', '-r200', out_png);
end
end
