function c3 = tf3d_from_design(row, p, opts)
%TF3D_FROM_DESIGN 3D coil shape and field, downstream of a chosen WP_TF design.
%
%   c3 = TF3D_FROM_DESIGN(row, p) takes one chosen design point (row, p as
%   for PLOT_WP_SECTION / WP_MECH_SURROGATE) and:
%     1. locates the current centroid of the winding pack on the inner leg
%        and mirrors it to the outer leg (fixing review point #2 of
%        coil3d/legacy: the original built the shape on the case face,
%        not the current centroid);
%     2. computes the bending-free shape of the coil's current centroid
%        line (TF_BENDING_FREE_SHAPE), analytic or iterated under the
%        real 3D field of all n_TF coils;
%     3. fits it with three tangent arcs (TF_THREE_ARC_FIT), the shape a
%        real coil is actually wound to;
%     4. builds a 3D filament model of every MADE turn on coil 1 (and one
%        filament per layer for the other n_TF-1 coils) and evaluates the
%        field with exact segment Biot-Savart (no ad hoc regularization
%        radius, review point #7);
%     5. checks the result against MADE's own 2D quantities: Ampere's law,
%        the field ripple, the 2D discrete field at the inner-leg midplane,
%        and the hoop tension T_bf;
%     6. optionally exports CAD offset curves (EXPORT_TF3D_SHAPE) and
%        plots (POSTPROCESS/PLOT_TF3D.M).
%
%   opts overrides p.tf3d_* (see input/WP_TF_input_template.xlsx, "3D
%   coil" category); all have defaults so this runs with an older input
%   file too.
%
%   c3 fields: geo (centroid radii, currents), shape (half D-shape, mode
%   0/1), arcs (three-arc fit), field (evaluation points and |B| on coil
%   1's turns), checks (struct of the comparisons in step 5).

if nargin < 3, opts = struct(); end
g = compute_operating_params(p);

%% 1. Geometry: current centroid and its mirror on the outer leg
field2d = compute_discrete_field_profile(row, p);   % also gives the validated 2D peak per turn
r_c = mean(field2d.y);                              % current centroid radius, inner leg (equal current per turn)
d_c = row.Ri_ - r_c;                                % centroid set-back from the case's plasma-facing face

geo.r_c = r_c; geo.d_c = d_c;
geo.r1 = g.RTFi - d_c;
geo.r2 = g.RTFo + d_c;
geo.n_TF = p.n_TF;
geo.n_turns_coil = sum(row.n_turns(1:row.n_layers));
geo.I_coil = row.Iop*geo.n_turns_coil;               % real total current per coil
geo.NI = g.NI;                                       % MADE's own total ampere-turns, for comparison

% Case envelope (for the CAD export), same recursion as everywhere else in
% MADE (search/scan_wp_designs.m, physics/wp_mech_surrogate.m surr_geometry)
Y1 = row.Ri_ - p.dr_plasma_side;                     % WP inner (plasma-side) radius
Y2 = Y1 - p.GoundIns - sum(row.Cond_h(1:row.n_layers)) - (row.n_layers-1)*p.INS_grades - p.GoundIns;
geo.WP_i = Y1; geo.WP_e = Y2;                        % inner leg: plasma-side / back of WP
geo.case_i = row.Ri_; geo.case_e = row.Rk_;          % inner leg: plasma-side face / nose inner radius
geo.nose_il = Y2 - row.Rk_;
geo.nose_ol_factor = get_or(p, opts, 'tf3d_nose_ol_factor', 0.5);
geo.nose_ol = geo.nose_ol_factor*geo.nose_il;

%% 2. Bending-free shape of the current centroid line
shape_opts.mode = get_or(p, opts, 'tf3d_shape_mode', 1);
shape_opts.n_seg = get_or(p, opts, 'tf3d_n_seg', 240);
shape_opts.tol = get_or(p, opts, 'tf3d_shape_tol', 1e-4);
shape_opts.max_iter = get_or(p, opts, 'tf3d_max_iter', 20);
shape = tf_bending_free_shape(geo.r1, geo.r2, geo.n_TF, geo.I_coil, shape_opts);

%% 3. Three-arc manufacturable fit
arcs = tf_three_arc_fit(shape);

%% 4. Filament model and field on coil 1's turns
use_arcs = get_or(p, opts, 'tf3d_use_arcs', 1);
n_fil_layer = get_or(p, opts, 'tf3d_n_fil_layer', 1);
n_seg_path = get_or(p, opts, 'tf3d_n_seg', 240);
if use_arcs
    path_rz = arc_path_points(arcs, n_seg_path);
else
    path_rz = [shape.r, shape.z];
end
field = tf3d_coil_filaments(row, p, geo, path_rz, n_fil_layer, field2d);

%% 5. Checks against MADE's own 2D/analytic quantities
Mu_0 = 4e-7*pi;
checks.ampere_ratio = field.B_phi_R0/(Mu_0*geo.n_TF*geo.I_coil/(2*pi*p.R0));
checks.ripple = field.ripple;
checks.ripple_target = p.ripple;
checks.T_axial = field.T_axial;
checks.T_bf_2d = 0.5*g.k_bf*geo.n_TF*(geo.I_coil)^2*Mu_0/(2*pi);
checks.T_ratio = checks.T_axial/checks.T_bf_2d;
checks.peak3d_mid = field.peak_mid;
checks.peak2d_mid = field.peak2d_mid;
checks.peak_ratio_mid = field.peak_mid/field.peak2d_mid;
checks.shape_mode_used = shape.mode;
checks.arc_fit_rms = arcs.rms;
checks.arc_fit_max_dev = arcs.max_dev;
checks.summary = sprintf(['Ampere ratio %.4f | ripple %.4f (target %.4f) | T_axial/T_bf(2D) %.4f | ' ...
    'peak 3D/2D at inner-leg midplane %.4f | arc fit RMS %.1f mm, max %.1f mm'], ...
    checks.ampere_ratio, checks.ripple, checks.ripple_target, checks.T_ratio, checks.peak_ratio_mid, ...
    1e3*arcs.rms, 1e3*arcs.max_dev);

c3.geo = geo; c3.shape = shape; c3.arcs = arcs; c3.field = field; c3.checks = checks;
c3.path_rz = path_rz;

%% 6. Optional CAD export
if get_or(p, opts, 'tf3d_export', 1)
    tag = get_or(p, opts, 'tf3d_tag', 'design');
    out_dir = get_or(p, opts, 'tf3d_out_dir', pwd);
    c3.export_files = export_tf3d_shape(c3, out_dir, tag);
end
end

function v = get_or(p, opts, name, default)
if isfield(opts, name) && ~isempty(opts.(name))
    v = opts.(name);
elseif isfield(p, name) && ~isempty(p.(name))
    v = p.(name);
else
    v = default;
end
end

function rz = arc_path_points(arcs, n_seg)
% Sample the three-arc fit into a single half-shape polyline (r,z), public
% convention: index 1 = A (inner leg top), index end = B (outer midplane).
n = max(20, round(n_seg/3));
P1 = sample_arc(arcs.C1, arcs.R1, arcs.ang1, n);
P2 = sample_arc(arcs.C2, arcs.R2, arcs.ang2, n);
P3 = sample_arc(arcs.C3, arcs.R3, arcs.ang3, n);
rz = [P1; P2(2:end,:); P3(2:end,:)];
end

function P = sample_arc(C, R, ang, n)
a0 = ang(1); a1 = ang(2);
d = a1 - a0;
if d > pi, d = d - 2*pi; elseif d < -pi, d = d + 2*pi; end
th = a0 + d*linspace(0, 1, n)';
P = [C(1) + R*cos(th), C(2) + R*sin(th)];
end
