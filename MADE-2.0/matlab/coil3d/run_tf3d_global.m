function G = run_tf3d_global(row, p, opts)
%RUN_TF3D_GLOBAL 3D beam + shell model of the whole TF system, MATLAB only.
%
%   G = RUN_TF3D_GLOBAL(row, p, opts) runs the global model of
%   TF3D_GLOBAL_MODEL / TF3D_GLOBAL_SOLVE on a MADE design point, without
%   ANSYS:
%     1. centreline of coil 1 from the 3D module (TF3D_FROM_DESIGN:
%        bending-free shape, optionally iterated on the discrete coil set);
%     2. nodal Lorentz forces on it from TF3D_CENTRELINE_LOADS (all n_TF
%        coils, rectangular WP section; PF/CS/plasma rings in opts.pf);
%     3. beams along the centreline of every coil + shells of the inner-leg
%        vault (and of the outer intercoil structure if opts.ois_R is given),
%        gravity support of every coil; linear static solution.
%
%   row, p: a design point and the machine parameters, as returned by
%   LOAD_DESIGN_POINT or chosen in MAIN_WP_TF_DESIGN.
%   opts (all optional):
%     W, t_wall   hollow square beam section of the case [m]
%                 (default: W = case width 2*Ri_*tan(pi/n_TF), t_wall 0.05;
%                 PROVISIONAL until the section is derived from the 2D design)
%     t_shell     shell thickness [m] (default 0.14, as STR_360)
%     ois_R       [Rmin Rmax] of the outer intercoil structure [m]
%                 (default [] = no OIS: machine-specific, give it explicitly)
%     pf          K x 5 [rc zc dr dz I] PF/CS/plasma rings (out-of-plane loads)
%     n_seg       segments of the centreline (default from TF3D_DEFAULTS)
%     tf3d        struct of options passed to TF3D_FROM_DESIGN
%     plot        true (default) to plot forces/moments along the coil
%
%   G: model m, solution res, loads L, centreline P, summary s (MN, MN m, mm).
%
%   Example (from MADE-2.0/matlab):
%     [row, p] = load_design_point('<tag>_results_<date>.xlsx', idx);
%     G = run_tf3d_global(row, p, struct('ois_R', [3.5 5]));

if nargin < 3, opts = struct(); end
tf3d_opts = getf(opts, 'tf3d', struct());
if isfield(opts, 'n_seg'), tf3d_opts.n_seg = opts.n_seg; end
if ~isstruct(row), row = table2struct(row); end

% 1. centreline of coil 1 (plane x = 0, y major radius, z vertical)
t3 = tf3d_from_design(row, p, tf3d_opts);
rz = t3.loop.rz;
if norm(rz(end,:) - rz(1,:)) < 1e-9, rz(end,:) = []; end   % closed loop: drop the repeated point
P = [zeros(size(rz,1),1), rz(:,1), rz(:,2)];
d = t3.design;

% 2. loads
L = tf3d_centreline_loads(P, d.wp_depth, d.wp_width, d.NI, p.n_TF, getf(opts, 'pf', zeros(0,5)));

% 3. global model
Ri = row.Ri_;
mo = struct('n_TF', p.n_TF, ...
    'W', getf(opts, 'W', 2*Ri*tan(pi/p.n_TF)), ...
    't_wall', getf(opts, 't_wall', 0.05), ...
    't_shell', getf(opts, 't_shell', 0.14));
ois_R = getf(opts, 'ois_R', []);
if isempty(ois_R)
    mo.with_ois = false;
else
    mo.with_ois = true; mo.ois_R = ois_R;
end
mo.gs_R = getf(opts, 'gs_R', Ri + 0.5*(max(P(:,2)) - Ri));   % gravity support search radius
m = tf3d_global_model(P, L.F, mo);
t0 = tic; res = tf3d_global_solve(m); ts = toc(t0);

% summary on coil 1 (beam forces: [N Vy Vz T My Mz] at node I)
M = size(P,1); bf = res.beam_f(1:M, :); U = res.U(1:M, :);
s = struct();
s.N_MN = [min(bf(:,1)) max(bf(:,1))]/1e6;
s.Mz_inplane_MNm = [min(bf(:,6)) max(bf(:,6))]/1e6;
s.My_outplane_MNm = [min(bf(:,5)) max(bf(:,5))]/1e6;
s.T_torsion_MNm = max(abs(bf(:,4)))/1e6;
s.U_mm = max(abs(U(:,1:3)), [], 1)*1e3;      % |ux| toroidal, |uy| radial, |uz| vertical
s.solve_s = ts; s.ndof = res.info.ndof; s.residual = res.info.residual;

fprintf('\nTF 3D global model: %d coils x %d nodes, %d beams, %d shells, %d DOF, %.1f s (residual %.1e)\n', ...
    p.n_TF, M, size(m.beam,1), size(m.shell,1), s.ndof, ts, s.residual);
fprintf('  section W %.3f m, wall %.3f m, shell %.3f m, OIS %s\n', mo.W, mo.t_wall, mo.t_shell, mat2str(ois_R));
fprintf('  axial force N        %8.2f .. %8.2f MN\n', s.N_MN);
fprintf('  in-plane moment      %8.2f .. %8.2f MN m\n', s.Mz_inplane_MNm);
fprintf('  out-of-plane moment  %8.2f .. %8.2f MN m\n', s.My_outplane_MNm);
fprintf('  max torsion          %8.2f MN m\n', s.T_torsion_MNm);
fprintf('  max |u| tor/rad/vert %6.2f / %6.2f / %6.2f mm\n', s.U_mm);

if getf(opts, 'plot', true)
    sc = [0; cumsum(sqrt(sum(diff(P).^2, 2)))]; sc = sc/sc(end);
    figure('Name', 'TF 3D global model');
    subplot(2,1,1); plot(sc, bf(:,1)/1e6); grid on; ylabel('N [MN]');
    title('Coil 1: axial force and moments along the centreline');
    subplot(2,1,2); plot(sc, bf(:,6)/1e6, sc, bf(:,5)/1e6, sc, bf(:,4)/1e6); grid on;
    legend('M in-plane', 'M out-of-plane', 'torsion'); xlabel('s / s_{max}'); ylabel('[MN m]');
end

G = struct('model', m, 'res', res, 'loads', L, 'P', P, 'summary', s, 'tf3d', t3, 'options', mo);
end

function v = getf(s, f, dflt)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = dflt; end
end
