function gm = tf3d_global_from_design(row, p, opts)
%TF3D_GLOBAL_FROM_DESIGN 3D beam + shell model of the whole TF system built from a MADE design point.
%
%   gm = TF3D_GLOBAL_FROM_DESIGN(row, p) builds and solves the global model
%   of TF3D_GLOBAL_MODEL / TF3D_GLOBAL_SOLVE (the MATLAB counterpart of the
%   ANSYS model STR_360, verified against it: docs/MODELLO_GLOBALE_3D.md)
%   with everything taken from the design point instead of the fixed
%   STR_360 data:
%
%   - centreline: the current centroid line of the WP (TF3D_DESIGN_GEOMETRY),
%     bending-free D from r1 (inner leg) to r2 (outer leg), straight inner
%     leg, n_seg segments (TF3D_RESAMPLE_LOOP); or the loop of a 3D step
%     already run (opts.tf3d = output of TF3D_FROM_DESIGN);
%   - loads: Lorentz forces at the centreline nodes from the 3D field of all
%     n_TF coils (TF3D_CENTRELINE_LOADS, as the ANSYS MAG_360 + its
%     post-processing), WP section n_gauss x n_gauss filaments; optional
%     PF/CS/plasma rings opts.pf (K x 5 [rc zc dr dz I]);
%   - beam section: the inner-leg section of the design, case + jackets
%     (TF3D_SECTION_FROM_DESIGN); vault shells as thick as the steel that
%     crosses the coil mid plane (nose + plate + jacket walls); OIS shells
%     t_ois (STR_360: 0.14 m);
%   - OIS zone and gravity support at the same fractions of the radial span
%     of the centreline as STR_360 (R 3.5-5 m and R <= 3 m on a span 1.082
%     - 5.207 m).
%
%   COUPLING WITH THE 2D SIZING. The 2D sizing (scan, layered model, case
%   and jacket surrogates, 2D FE) loads the inner leg with the
%   bending-free vertical force T_bf = 0.5 k_bf n_TF (N I)^2 mu0/(2 pi),
%   half of the vertical Lorentz force of the upper half of a coil. The
%   global model gives the force that the inner leg really carries:
%   beam axial force + vertical membrane force of the vault shells over
%   one chord, F_inner(z). Their ratio is the axial load factor
%
%     k_axial = max over the straight inner leg of F_inner / T_bf
%
%   (axial_mode = 1: F_inner + |M_inplane| c_r A/I_z, the extreme fibre).
%   With p.axial_load_factor = k_axial every sizing step uses k_axial T_bf
%   (AXIAL_LOAD_FACTOR); COUPLE_GLOBAL_SIZING iterates sizing and global
%   model to convergence.
%
%   opts (or Excel p.g3d_<name>): n_seg (160), n_gauss (4), t_ois (0.14 m),
%   ois_frac ([0.586 0.950]), gs_frac (0.465), axial_mode (0), pf, tf3d,
%   verbose (true);
%   ois ([] = the ois_frac window): K x 3 [s_start s_end t] OIS panels,
%   position as fractions of the centreline (0 = inner mid plane, ~0.5 =
%   outer mid plane) and thickness [m]; every OIS layout is checked by
%   TF3D_OIS_ZONES before the run (overlaps, vault, free width between the
%   cases, clearance ois_clearance (0.1 m) from the PF rings), and shown
%   and confirmed if ois_plot / ois_confirm are true.
%
%   gm: k_axial, k_axial_mid (at the mid plane), T_bf (k = 1), F_inner,
%   N_beam, N_vault (per inner-leg element), z_inner, M_inplane (coil 1,
%   per element), sigma (struct of stresses), Fz_half (vertical Lorentz
%   force of the upper half of coil 1), section, model, res, loads, P,
%   options, elapsed.

if nargin < 3, opts = struct(); end
t0 = tic;
o = struct('n_seg', 160, 'n_gauss', 4, 't_ois', 0.14, 'ois_frac', [0.586 0.950], 'gs_frac', 0.465, ...
    'axial_mode', 0, 'pf', zeros(0, 5), 'tf3d', [], 'verbose', true, ...
    'ois', [], 'ois_plot', false, 'ois_confirm', false, 'ois_clearance', 0.1);
fn = fieldnames(o);
for i = 1:numel(fn)
    key = ['g3d_' fn{i}];
    if isfield(p, key) && ~isempty(p.(key)), o.(fn{i}) = p.(key); end
end
fn = fieldnames(opts);
for i = 1:numel(fn)
    if ~isfield(o, fn{i}), error('tf3d_global_from_design:option', 'Unknown option: %s', fn{i}); end
    o.(fn{i}) = opts.(fn{i});
end
if ~isstruct(row), row = table2struct(row); end
if ~any(o.axial_mode == [0 1]), error('tf3d_global_from_design:option', 'axial_mode must be 0 or 1.'); end

% centreline of coil 1 (plane x = 0, y = major radius, z vertical)
d = tf3d_design_geometry(row, p);
if ~isempty(o.tf3d)
    loop = o.tf3d.loop;
else
    shape = tf_bending_free_shape(d.r1, d.r2, max(401, 2*o.n_seg + 1));
    loop = tf3d_resample_loop(shape.upper, o.n_seg);
end
rz = loop.rz(1:end-1, :);                    % closed loop: last point = first
M = size(rz, 1);
P = [zeros(M, 1), rz];

% loads: field of all coils at the centreline, nodal Lorentz forces
nl = row.n_layers;
A_WP = sum(row.n_turns(1:nl).*row.Cond_w(1:nl).*row.Cond_h(1:nl));
WPH = d.wp_depth; WPW = A_WP/WPH;           % equivalent rectangle of the stepped WP

% section and model
sec = tf3d_section_from_design(row, p);
R1 = rz(:, 1); Rmin = min(R1); Rmax = max(R1);
mo = struct('n_TF', p.n_TF, 'E', sec.E, 'nu', sec.nu, 'sec', sec, 't_vault', sec.t_vault, ...
    't_ois', o.t_ois, 'ois_R', Rmin + o.ois_frac*(Rmax - Rmin), 'gs_R', Rmin + o.gs_frac*(Rmax - Rmin), ...
    'vault_tol', 1e-6*Rmax);
% OIS panels: o.ois = [s_start s_end t; ...] (fractions of the centreline,
% 0 = inner mid plane, ~0.5 = outer mid plane), or the default R window
% ois_frac (upper and lower part, thickness t_ois); checked before the run
% (overlaps, vault, free width between the cases, PF clearance), and shown
% / confirmed with ois_plot / ois_confirm (TF3D_OIS_ZONES)
inn0 = find(abs(R1 - Rmin) < 1e-6*Rmax);
ois_spec = o.ois;
if isempty(ois_spec)
    s_c = [0; cumsum(sqrt(sum(diff(P).^2, 2)))]; s_c = s_c/(s_c(end) + norm(P(1,:) - P(end,:)));
    inO = R1 >= mo.ois_R(1) & R1 <= mo.ois_R(2);
    for half = [1 -1]
        k = find(inO & half*rz(:, 2) >= 0);
        if numel(k) >= 2, ois_spec = [ois_spec; s_c(min(k)) s_c(max(k)) o.t_ois]; end %#ok<AGROW>
    end
end
mo.ois_zones = tf3d_ois_zones(P, ois_spec, struct('n_TF', p.n_TF, 'W', 2*row.Ri_*tan(pi/p.n_TF), ...
    'vault_idx', inn0, 'pf', o.pf, 'clearance', o.ois_clearance, 'plot', o.ois_plot, 'confirm', o.ois_confirm));
% loads (after the OIS check: the slow step)
lo = tf3d_centreline_loads(P, WPH, WPW, d.NI, p.n_TF, o.pf, struct('n_gauss', o.n_gauss));
m = tf3d_global_model(P, lo.F, mo);
res = tf3d_global_solve(m);

% inner leg of coil 1: beam axial force + vault shells over one chord
inn = abs(R1 - Rmin) < 1e-6*Rmax;
nxt = [2:M 1];
ie = find(inn & inn(nxt));                    % elements j -> j+1 on the straight leg
z_e = (rz(ie, 2) + rz(nxt(ie), 2))/2;
N_beam = res.beam_f(ie, 1);
sc = m.shell_coil; sr = m.shell_region;
vs = find(sc == 1 & sr == 1);                 % vault shells between coil 1 and coil 2
% shell k of coil 1 spans nodes [j, j+1] of coil 1 (element j)
j_of = mod(m.shell(vs, 1) - 1, M) + 1;
N22 = nan(numel(ie), 1); N11 = N22; t_v = sec.t_vault;
smax_v = N22;
for q = 1:numel(ie)
    k = vs(j_of == ie(q));
    if isempty(k), continue, end
    r = res.shell_r(k(1), :);
    N11(q) = r(1); N22(q) = r(2);
    s_top = [r(1) r(2)]/t_v + 6*[r(4) r(5)]/t_v^2; s_bot = [r(1) r(2)]/t_v - 6*[r(4) r(5)]/t_v^2;
    smax_v(q) = max([tresca2(s_top, r(3)/t_v + 6*r(6)/t_v^2), tresca2(s_bot, r(3)/t_v - 6*r(6)/t_v^2)]);
end
chord = 2*Rmin*sin(pi/p.n_TF);
N_vault = N22*chord;
N22(isnan(N22)) = 0; N_vault(isnan(N_vault)) = 0;
F_inner = N_beam + N_vault;
Mz_in = max(abs(res.beam_f(ie, [6 12])), [], 2);
F_fibre = F_inner + Mz_in*sec.c_r*sec.A/sec.Iz;

g = compute_operating_params(p);
T_bf = 0.5*g.k_bf*p.n_TF*d.NI^2*(4e-7*pi)/(2*pi);       % 2D sizing, k = 1
if o.axial_mode == 1, Fk = F_fibre; else, Fk = F_inner; end
[k_max, iq] = max(Fk/T_bf);
[~, imid] = min(abs(z_e));
k_mid = Fk(imid)/T_bf;

% stresses of the beam (coil 1) with the design section
Bf = res.beam_f(1:M, :);
sN = Bf(:, 1)/sec.A;
sMz = max(abs(Bf(:, [6 12])), [], 2)*sec.c_r/sec.Iz;
sMy = max(abs(Bf(:, [5 11])), [], 2)*sec.c_t/sec.Iy;
tau_T = max(abs(Bf(:, [4 10])), [], 2)/(2*sec.A_m*sec.t_case_mean);
sigma = struct('axial', sN, 'bend_inplane', sMz, 'bend_outplane', sMy, 'torsion_shear', tau_T, ...
    'fibre_max', max(sN + sMz + sMy), 'inner_axial_max', max(F_inner)/sec.A, ...
    'inner_bend_max', max(Mz_in)*sec.c_r/sec.Iz, 'vault_tresca_max', max(smax_v), ...
    'vault_N11_range', [min(N11) max(N11)], 'vault_N22_range', [min(N22) max(N22)]);

% vertical Lorentz force of the upper half of coil 1 (nodes at z = 0 count half)
zz = rz(:, 2); wz = (zz > 1e-9) + 0.5*(abs(zz) <= 1e-9);
Fz_half = sum(wz.*lo.F(:, 3));

gm = struct('k_axial', k_max, 'k_axial_mid', k_mid, 'axial_mode', o.axial_mode, 'T_bf', T_bf, ...
    'F_inner', F_inner, 'F_fibre', F_fibre, 'N_beam', N_beam, 'N_vault', N_vault, 'z_inner', z_e, ...
    'z_kmax', z_e(iq), 'M_inplane_inner', Mz_in, 'vault_N11', N11, 'vault_N22', N22, 'vault_tresca', smax_v, ...
    'sigma', sigma, 'Fz_half', Fz_half, 'section', sec, 'model', m, 'res', res, 'loads', lo, 'P', P, ...
    'inner_elements', ie, 'options', o, 'design', d, 'elapsed', toc(t0));
if o.verbose
    print_summary(gm, p);
end
end

function s = tresca2(sp, txy)
% Tresca intensity of a plane stress state (s1, s2, t12), third stress 0
c = (sp(1) + sp(2))/2; r = hypot((sp(1) - sp(2))/2, txy);
s = max([2*r, abs(c + r), abs(c - r)]);
end

function print_summary(gm, p)
s = gm.section; sg = gm.sigma;
fprintf(['\nTF global model (beam + shell, all %d coils): %d nodes per coil, %d beams, %d shells, ' ...
    '%d DOF, %.1f s\n'], p.n_TF, size(gm.P, 1), size(gm.model.beam, 1), size(gm.model.shell, 1), ...
    gm.res.info.ndof, gm.elapsed);
fprintf(['  section of the design: steel A %.4f m^2 (case %.4f, jackets %.4f), I in-plane %.3g m^4, ' ...
    'out-of-plane %.3g m^4, J %.3g m^4, vault shell %.0f mm\n'], s.A, s.A_case, s.A_jacket, s.Iz, s.Iy, s.J, ...
    1e3*s.t_vault);
fprintf('  vertical Lorentz force of the upper half: %.2f MN = 2 x %.2f MN; 2D sizing T_bf = %.2f MN per leg\n', ...
    gm.Fz_half/1e6, gm.Fz_half/2e6, gm.T_bf/1e6);
fprintf(['  inner leg: axial force %.2f ... %.2f MN (beam %.2f ... %.2f + vault shells %.2f ... %.2f), ' ...
    'in-plane moment up to %.2f MN m\n'], min(gm.F_inner)/1e6, max(gm.F_inner)/1e6, min(gm.N_beam)/1e6, ...
    max(gm.N_beam)/1e6, min(gm.N_vault)/1e6, max(gm.N_vault)/1e6, max(gm.M_inplane_inner)/1e6);
fprintf(['  axial load factor k = F_inner/T_bf: %.3f (max, z = %.2f m), %.3f at the mid plane%s\n'], ...
    gm.k_axial, gm.z_kmax, gm.k_axial_mid, ifelse(gm.axial_mode == 1, ' (extreme fibre)', ''));
fprintf(['  inner leg stresses: axial %.0f MPa, in-plane bending %.0f MPa (extreme fibre); coil max fibre ' ...
    '%.0f MPa; vault: hoop N11 %.1f ... %.1f MN/m, Tresca %.0f MPa\n'], sg.inner_axial_max/1e6, ...
    sg.inner_bend_max/1e6, sg.fibre_max/1e6, sg.vault_N11_range/1e6, sg.vault_tresca_max/1e6);
end

function v = ifelse(c, a, b)
if c, v = a; else, v = b; end
end
