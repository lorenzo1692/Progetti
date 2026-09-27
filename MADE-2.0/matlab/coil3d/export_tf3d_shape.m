function files = export_tf3d_shape(c3, out_dir, tag)
%EXPORT_TF3D_SHAPE Write CAD offset curves and shape details for a 3D coil.
%
%   files = EXPORT_TF3D_SHAPE(c3, out_dir, tag) writes, from the coil's
%   current-centroid path (c3.path_rz, half-shape, public convention),
%   offset curves at constant perpendicular distance from it: the WP
%   plasma-side and back faces (constant thickness along the whole coil,
%   dr_plasma_side and WP_h from the chosen design point) and the case
%   plasma-facing and nose surfaces (nose thickness tapering linearly from
%   the inner-leg value geo.nose_il to the outer-leg value geo.nose_ol -
%   TO BE CONFIRMED, see docs/HANDOFF.md #3 - since only the inner leg is
%   sized by the WP_TF scan). Same file layout as coil3d/legacy (r,z pairs
%   in mm, one column pair per curve), plus a plain-text arc summary.
%
%   This replaces coil3d/legacy's fixed WP_w x WP_h=1.0x0.5 m ANSYS source
%   (review point #3) with the actual chosen design's dimensions, and
%   removes the dependency on a specific dated output filename.

if ~isfolder(out_dir), mkdir(out_dir); end
rz = c3.path_rz;                                     % public convention half-shape (centroid)
geo = c3.geo;

[nr, nz] = local_normals(rz);
n = size(rz, 1);
d_nose = linspace(geo.nose_il, geo.nose_ol, n)';      % linear taper, inner leg -> outer leg

curve_WP_i    = offset_curve(rz, nr, nz,  geo.WP_i - geo.r_c);
curve_WP_e    = offset_curve(rz, nr, nz,  geo.WP_e - geo.r_c);
curve_case_i  = offset_curve(rz, nr, nz,  geo.case_i - geo.r_c);
curve_case_e  = [rz(:,1) + nr.*(geo.WP_e - geo.r_c - d_nose), rz(:,2) + nz.*(geo.WP_e - geo.r_c - d_nose)];
curve_CL      = rz;

mirror = @(C) [ [C; flipud([C(1:end-1,1), -C(1:end-1,2)])] ];  %#ok<NBRAK> full closed curve, up-down mirrored

mm = 1e3;
shape_file = fullfile(out_dir, sprintf('%s_TF3D_shape_%s.txt', tag, datestr(now, 'yyyymmdd_HHMMSS'))); %#ok<TNOW1,DATST>
fid = fopen(shape_file, 'w');
fprintf(fid, 'r_case_CL z_case_CL r_case_i z_case_i r_case_e z_case_e r_WP_i z_WP_i r_WP_e z_WPe\n');
data = [mirror(curve_CL), mirror(curve_case_i), mirror(curve_case_e), mirror(curve_WP_i), mirror(curve_WP_e)]*mm;
fprintf(fid, '%f %f %f %f %f %f %f %f %f %f\n', data');
fclose(fid);

detail_file = fullfile(out_dir, sprintf('%s_TF3D_arcs_%s.txt', tag, datestr(now, 'yyyymmdd_HHMMSS'))); %#ok<TNOW1,DATST>
fid = fopen(detail_file, 'w');
fprintf(fid, 'TF 3D coil - three-arc fit (coil3d/tf_three_arc_fit.m)\n');
fprintf(fid, 'Shape mode used: %d (0=analytic bending-free, 1=iterated under the real n_TF-coil field)\n', c3.shape.mode);
if isfield(c3.shape, 'iter'), fprintf(fid, 'Shape iterations: %d (err_r=%.2e, err_z=%.2e)\n', c3.shape.iter, c3.shape.err_r, c3.shape.err_z); end
fprintf(fid, '\nC1: %.6f %.6f\nR1: %.6f\n', c3.arcs.C1, c3.arcs.R1);
fprintf(fid, '\nC2: %.6f %.6f\nR2: %.6f\n', c3.arcs.C2, c3.arcs.R2);
fprintf(fid, '\nC3: %.6f %.6f\nR3: %.6f\n', c3.arcs.C3, c3.arcs.R3);
fprintf(fid, '\nFit quality: RMS %.2f mm, max %.2f mm\n', 1e3*c3.arcs.rms, 1e3*c3.arcs.max_dev);
fprintf(fid, '\n%s\n', c3.checks.summary);
fclose(fid);

files = {shape_file, detail_file};
end

function [nr, nz] = local_normals(rz)
% Outward unit normal (increasing r on the inner leg) at every point of
% the half-shape, from a central-difference tangent.
n = size(rz, 1);
t = zeros(n, 2);
t(2:end-1,:) = rz(3:end,:) - rz(1:end-2,:);
t(1,:) = rz(2,:) - rz(1,:);
t(end,:) = rz(end,:) - rz(end-1,:);
len = sqrt(sum(t.^2, 2)); t = t./len;
nr = -t(:,2); nz = t(:,1);                            % rotate tangent by +90 deg
if nr(1) < 0, nr = -nr; nz = -nz; end                 % fix sign: outward = increasing r at the inner leg
end

function C = offset_curve(rz, nr, nz, d)
C = [rz(:,1) + nr*d, rz(:,2) + nz*d];
end
