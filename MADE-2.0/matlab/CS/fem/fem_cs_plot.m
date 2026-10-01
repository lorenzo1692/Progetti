function fem_cs_plot(res, c)
%FEM_CS_PLOT Plot the hoop stress of load case c on the stack and the axial force through the plates.
%
%   FEM_CS_PLOT(res, c) draws the winding-pack elements of the stack colored
%   by the jacket hoop stress [MPa] and, beside it, the compression through
%   every spacer plate [MN] for load case c of the result of FEM_CS_VERIFY.

if nargin < 2, c = 1; end
st = res.stack; cs = res.cases(c);
ew = find(st.matid == 1);
figure('units', 'normalized', 'outerposition', [0.1 0.1 0.8 0.85]);
subplot(1, 2, 1);
mesh = st.mesh;
patch('Faces', mesh.elems, 'Vertices', mesh.nodes, 'FaceColor', [0.85 0.85 0.85], 'EdgeColor', 'none'); hold on
patch('Faces', mesh.elems(ew, :), 'Vertices', mesh.nodes, 'FaceVertexCData', cs.s_th*1e-6, ...
    'FaceColor', 'flat', 'EdgeColor', 'none');
colormap(jet); cb = colorbar; cb.Label.String = 'jacket hoop stress [MPa]';
axis equal; xlabel('r [m]'); ylabel('z [m]');
title(sprintf('CS stack - %s', cs.name), 'Interpreter', 'none');
subplot(1, 2, 2);
ip = find(st.aband == 3);
zp = 0.5*(st.zv(ip) + st.zv(ip+1));
plot(cs.F_plate_MN, zp, 'o-', 'LineWidth', 1.5); grid on
xlabel('compression through the plate [MN]'); ylabel('z [m]');
title('Axial force along the stack');
end
