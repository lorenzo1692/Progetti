function shape = tf_bending_free_shape(r1, r2, n_TF, I_coil, opts)
%TF_BENDING_FREE_SHAPE Bending-free ("Princeton D") shape of a TF coil.
%
%   shape = TF_BENDING_FREE_SHAPE(r1, r2, n_TF, I_coil, opts) returns the
%   in-plane (r,z) shape of one TF coil's CURRENT CENTROID line that is in
%   pure tension (no bending moment) under the toroidal field, from the
%   inner-leg centroid radius r1 to the outer-leg centroid radius r2.
%
%   This fixes review point #2 of the original scripts (coil3d/legacy):
%   they built the bending-free shape on the case's plasma-facing radii
%   (Ri, Re) instead of on the current centroid, which is the only line
%   the bending-free condition actually applies to.
%
%   opts.mode (from p.tf3d_shape_mode, default 1):
%     0 - analytic closed form (the classical Princeton-D / File shape),
%         assuming the idealized continuous-sheet field
%         B_eff(r) = mu0*n_TF*I_coil/(2*pi*r). Always converges, ignores
%         field ripple from the finite number of coils.
%     1 - iterate the shape under the REAL field of all n_TF coils: the
%         analyzed coil is an ideal filament of its own (self term
%         mu0*I_coil/(2*pi*r), the usual idealization for a coil's own
%         hoop-tension balance) plus the exact Biot-Savart field of the
%         other n_TF-1 coils (each the current shape estimate, rotated by
%         its toroidal angle), evaluated ON the coil's own candidate path.
%         This is what coil3d/legacy/bendingfree_opt.m did by calling
%         ANSYS every iteration (review point #1: not portable, needs a
%         Windows machine and a licensed solver); the same Newton
%         iteration is kept here (verified to reproduce the mode-0 closed
%         form when fed the same idealized field, see
%         VALIDATION/VALIDATE_TF3D.M), only the field evaluation is
%         replaced by BIOT_SAVART_SEGMENTS. Falls back to mode 0 with a
%         warning if it does not converge within opts.max_iter.
%
%         Calibration note: a point sitting exactly on the current sheet
%         formed by n_TF coils sees, from the OTHER n_TF-1 coils alone,
%         only HALF of the idealized "deep inside the bore" continuum
%         value mu0*n_TF*I/(2*pi*r) as n_TF grows (the field has a jump
%         discontinuity across the sheet, and its value exactly on the
%         sheet is the average of the two sides, here bore-side and
%         zero outside) - verified numerically in
%         VALIDATION/VALIDATE_TF3D.M. MADE's own 2D model (B_PHI_TF,
%         T_bf, ...) uses the non-averaged "deep inside" convention
%         throughout, so the discrete neighbor contribution here is
%         doubled to match it (self stays the plain ideal-filament term,
%         already 1/n_TF of the continuum value and negligible at
%         realistic n_TF): this makes mode 1 converge to mode 0 as
%         n_TF -> infinity by construction, while still capturing the
%         real ripple at the machine's actual n_TF.
%
%   opts.n_seg (from p.tf3d_n_seg, default 240): points describing the
%   half-shape; opts.tol (from p.tf3d_shape_tol, default 1e-4): relative
%   change in r and z between iterations below which mode 1 stops;
%   opts.max_iter (from p.tf3d_max_iter, default 20).
%
%   shape fields (public convention: index 1 = inner leg, index end =
%   outer leg, opposite of the internal Newton-solve convention below):
%     t        - parameter in [-pi/2, pi/2], t(1) = inner leg
%     r, z     - half D-shape, from the top of the inner straight leg
%                (r1, z(1) > 0) to the outer-leg midplane point (r2, 0)
%     mode     - 0 or 1, the mode actually used (1 may fall back to 0)
%     iter, err_r, err_z - convergence info (mode 1 only)

if nargin < 5, opts = struct(); end
mode = get_opt(opts, 'mode', 1);
N = get_opt(opts, 'n_seg', 240);

if mode == 0
    shape = to_public(newton_shape_ideal(r1, r2, n_TF, I_coil, N));
    shape.mode = 0;
    return
end

tol = get_opt(opts, 'tol', 1e-4);
max_iter = get_opt(opts, 'max_iter', 20);
Mu_0 = 4e-7*pi;

% Internal (Newton-solve) convention: t=-pi/2 -> r2 (outer, z=0),
% t=+pi/2 -> r1 (inner leg top); verified against the mode-0 closed form
% under the idealized field (see VALIDATION/VALIDATE_TF3D.M).
%
% Initial guess: the analytic (mode-0) shape, not a crude circular arc.
% Unlike the ideal 1/r field (a fixed function of r alone, for which even
% a circular-arc start converges cleanly), Beff here is Biot-Savart from
% the CANDIDATE SHAPE ITSELF, so a poor starting shape can converge to a
% spurious nearby fixed point of the self-referential system instead of
% the physical one. Starting from mode 0 (already exact in the large-n_TF
% limit, and verified to be a stable fixed point of this same iteration
% at n_TF=200 in VALIDATION/VALIDATE_TF3D.M) avoids that.
t = linspace(-pi/2, pi/2, N)';
[r, z] = newton_shape_ideal_rz(r1, r2, n_TF, I_coil, t);

converged = false; it = 0;
for it = 1:max_iter
    [S1, S2] = closed_path_segments(r, z);          % this coil's own closed path (source)
    P = [r, zeros(size(r)), z];                      % evaluation points, phi = 0 frame
    Beff = Mu_0*I_coil./(2*pi*r);                     % self, ideal-filament term
    for c = 2:n_TF
        a = (c-1)*2*pi/n_TF;
        Sc1 = rot_toroidal(S1, a); Sc2 = rot_toroidal(S2, a);
        Bc = biot_savart_segments(P, Sc1, Sc2, I_coil);
        Beff = Beff + 2*Bc(:,2);                       % toroidal component = y at phi=0; x2: sheet-average calibration above
    end
    [r_new, z_new] = newton_step(t, r, Beff);

    err_r = norm(r_new - r)/norm(r); err_z = norm(z_new - z)/max(norm(z), eps);
    r = r_new; z = z_new;
    if err_r < tol && err_z < tol, converged = true; break, end
end

if ~converged
    warning('tf_bending_free_shape:not_converged', ...
        ['3D shape iteration did not converge in %d steps (err_r=%.2e, err_z=%.2e): ' ...
         'falling back to the analytic bending-free shape (mode 0).'], max_iter, err_r, err_z);
    shape = to_public(newton_shape_ideal(r1, r2, n_TF, I_coil, N));
    shape.mode = 0; shape.iter = it; shape.err_r = err_r; shape.err_z = err_z;
    return
end
shape = to_public(struct('t', t, 'r', r, 'z', z));
shape.mode = 1; shape.iter = it; shape.err_r = err_r; shape.err_z = err_z;
end

function v = get_opt(opts, name, default)
if isfield(opts, name) && ~isempty(opts.(name)), v = opts.(name); else, v = default; end
end

function shp = newton_shape_ideal(r1, r2, n_TF, I_coil, N)
% Mode-0 analytic shape, computed by the SAME Newton solver as mode 1 but
% fed the idealized continuous-sheet field: this guarantees mode 0 and
% mode 1 agree exactly in the large-n_TF / no-ripple limit by construction
% (rather than by two independently-coded formulas that could drift).
t = linspace(-pi/2, pi/2, N)';
[r, z] = newton_shape_ideal_rz(r1, r2, n_TF, I_coil, t);
shp = struct('t', t, 'r', r, 'z', z);
end

function [r, z] = newton_shape_ideal_rz(r1, r2, n_TF, I_coil, t)
% r,z (internal convention) solving the bending-free equations under the
% idealized continuous-sheet field mu0*n_TF*I/(2*pi*r); a fixed function
% of r alone, so a circular-arc start converges reliably (unlike mode 1's
% shape-dependent Biot-Savart field, see the comment where this is called).
Mu_0 = 4e-7*pi;
r0 = (r2 + r1)/2 + (r2 - r1)/2*cos(t + pi/2);
Beff = Mu_0*n_TF*I_coil./(2*pi*r0);
[r, z] = newton_step(t, r0, Beff);
% one more pass to tighten (Beff depends only on r, already close after one step)
Beff = Mu_0*n_TF*I_coil./(2*pi*r);
[r, z] = newton_step(t, r, Beff);
end

function [r_new, z_new, k] = newton_step(t, r, Beff)
% One pass of the bending-free Newton update (r-solve then z-integration)
% given the field Beff(t) already evaluated at the current (r,z). Internal
% convention: t=-pi/2 -> outer leg (r2, z=0), t=+pi/2 -> inner leg (r1).
N = numel(t);
rr = flipud(r); Bs = flipud(Beff);                    % increasing-r order, for the r1->r table
intB = zeros(N, 1);
for i = 2:N
    intB(i) = intB(i-1) + 0.5*(Bs(i-1) + Bs(i))*(rr(i) - rr(i-1));
end
k = intB(end)/2;
r_new = r;
for i = 1:N
    x0 = r(i);
    for nit = 1:60
        F = interp1(rr, intB, x0, 'linear', 'extrap') + k*(sin(t(i)) - 1);
        dF = interp1(rr, Bs, x0, 'linear', 'extrap');
        xx = x0 - F/dF;
        conv = abs(xx - x0)/abs(x0) < 1e-12;
        x0 = xx;
        if conv, break, end
    end
    r_new(i) = x0;
end
z_new = zeros(N, 1);
for i = 2:N
    z_new(i) = z_new(i-1) - k*(sin(t(i-1))/Beff(i-1) + sin(t(i))/Beff(i))/2*(t(i) - t(i-1));
end
end

function pub = to_public(internal)
% Flip the internal (outer-first) convention to the public (inner-first,
% z(1) > 0) one documented at the top of this file.
pub.t = -flipud(internal.t);
pub.r = flipud(internal.r);
pub.z = flipud(internal.z);
end

function [S1, S2] = closed_path_segments(r, z)
% Full closed coil path (one period) built from the half-shape in the
% INTERNAL convention (r(1)=r2,z=0 .. r(end)=r1,z=z_top): upper half as
% given, mirrored lower half, and the straight inner leg at r=r(end)
% connecting the two ends.
ru = flipud(r); zu = flipud(z);                       % now r1 (top) .. r2 (z=0), increasing z along the path
rl = flipud(ru); zl = -flipud(zu);                    % mirror: r2 (z=0) .. r1 (z=-z_top)
rs = [ru; rl(2:end)];
zs = [zu; zl(2:end)];
rs = [rs; rs(1)]; zs = [zs; zs(1)];                   % close with the straight inner leg
x = rs; y = zeros(size(rs)); zc = zs;
S1 = [x(1:end-1), y(1:end-1), zc(1:end-1)];
S2 = [x(2:end),   y(2:end),   zc(2:end)];
end

function Sr = rot_toroidal(S, a)
% Rotate 3D points by angle a about the z axis (toroidal placement of
% another coil's identical shape).
c = cos(a); s = sin(a);
Sr = [S(:,1)*c - S(:,2)*s, S(:,1)*s + S(:,2)*c, S(:,3)];
end
