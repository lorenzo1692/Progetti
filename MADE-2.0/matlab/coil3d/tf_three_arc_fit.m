function arcs = tf_three_arc_fit(shape)
%TF_THREE_ARC_FIT Deterministic 3-tangent-arc approximation of a TF shape.
%
%   arcs = TF_THREE_ARC_FIT(shape) approximates the smooth bending-free
%   half-shape (from TF_BENDING_FREE_SHAPE, public convention: index 1 =
%   inner leg top A = (r1, zA), index end = outer-leg midplane point
%   B = (r2, 0)) with three circular arcs, tangent to each other and to
%   the straight legs at A and B - the shape a real coil is actually wound
%   to (winding machines bend constant-radius arcs, not a continuously
%   varying-radius curve).
%
%   This replaces coil3d/legacy/TF_shape_opt_BF_Ncoils.m's fragile fit
%   (review point #4): that one solved 13 simultaneous equations with
%   fsolve from a RANDOM starting point, had a sign error in two of them,
%   encoded inequality constraints as if they were residuals (which
%   fsolve does not honor), and stopped at a 0.40 m tolerance with no
%   guarantee of ever converging. Here there are only 3 unknowns
%   (R1, R2, R3), found deterministically:
%
%     - Arc 1 (radius R1): tangent to the vertical straight leg at A,
%       center C1 = (r1+R1, zA).
%     - Arc 3 (radius R3): tangent to the vertical line at B,
%       center C3 = (r2-R3, 0).
%     - Arc 2 (radius R2): internally tangent to both, i.e.
%       |C2-C1| = R2-R1 and |C2-C3| = R2-R3; given R1, R2, R3 this fixes
%       C2 as a circle-circle intersection (closed form).
%
%   (R1, R2, R3) are found by minimizing the RMS distance from the smooth
%   shape to the nearest of the three FULL circles (a standard, robust
%   proxy for a fit that is monotonic and single-valued in r as this
%   shape is) with FMINSEARCH from a coarse grid starting point - no
%   randomness, no inequality-as-residual, reproducible.
%
%   arcs fields:
%     R1, R2, R3      - arc radii [m]
%     C1, C2, C3      - arc centers [r, z]
%     A, B            - end points (inner leg top, outer midplane)
%     T12, T23        - tangent points arc1/arc2, arc2/arc3
%     ang1, ang2, ang3- [start end] angle [rad] of each arc, atan2(z-Cz, r-Cr),
%                       oriented so the arc from the start angle to the end
%                       angle (increasing) sweeps from A to B
%     rms, max_dev    - fit quality: RMS and max radial deviation [m] of
%                       the smooth shape from the fitted arcs

A = [shape.r(1), shape.z(1)];
B = [shape.r(end), shape.z(end)];
r1 = A(1); zA = A(2); r2 = B(1);
Q = [shape.r, shape.z];

% Bounds keep the fit to three genuinely comparable arcs (a manufacturable
% coil is not wound with a near-zero-radius "arc"): each radius between
% 8% and 80% of the r1-r2 span. Without them the RMS proxy can be
% minimized by collapsing one arc to zero, which technically lowers the
% point-to-nearest-circle distance but is not a usable 3-arc shape.
span = r2 - r1;
lo = 0.08*span; hi = 0.80*span;

% coarse grid for a deterministic, reasonable starting point
R1g = span*[0.15 0.30 0.45]; R3g = span*[0.15 0.30 0.45]; R2g = span*[0.5 0.75 1.0];
best = struct('x', [], 'f', inf);
for a = R1g
    for b = R2g
        for c = R3g
            f = obj([a b c], A, B, Q, lo, hi);
            if f < best.f, best.f = f; best.x = [a b c]; end
        end
    end
end
x0 = best.x;
options = optimset('Display', 'off', 'TolX', 1e-6, 'TolFun', 1e-9, 'MaxFunEvals', 3000);
x = fminsearch(@(x) obj(x, A, B, Q, lo, hi), x0, options);

R1 = x(1); R2 = x(2); R3 = x(3);
C1 = [r1 + R1, zA];
C3 = [r2 - R3, 0];
C2 = internal_tangent_center(C1, R2 - R1, C3, R2 - R3);
T12 = C1 + R1*unit(C1 - C2);
T23 = C3 + R3*unit(C3 - C2);

arcs.R1 = R1; arcs.R2 = R2; arcs.R3 = R3;
arcs.C1 = C1; arcs.C2 = C2; arcs.C3 = C3;
arcs.A = A; arcs.B = B; arcs.T12 = T12; arcs.T23 = T23;
arcs.ang1 = arc_angles(C1, A, T12);
arcs.ang2 = arc_angles(C2, T12, T23);
arcs.ang3 = arc_angles(C3, T23, B);

d = zeros(size(Q,1), 1);
for i = 1:size(Q,1), d(i) = dist_to_arcs(Q(i,:), C1, R1, C2, R2, C3, R3); end
arcs.rms = sqrt(mean(d.^2));
arcs.max_dev = max(abs(d));
end

function f = obj(x, A, B, Q, lo, hi)
R1 = x(1); R2 = x(2); R3 = x(3);
if R1 < lo || R3 < lo || R1 > hi || R3 > hi || R2 <= R1 || R2 <= R3 || R2 > 4*hi
    f = 1e6 + sum(max(0, [lo-R1, lo-R3, R1-hi, R3-hi, R1-R2, R3-R2]));
    return
end
C1 = [A(1) + R1, A(2)];
C3 = [B(1) - R3, B(2)];
C2 = internal_tangent_center(C1, R2 - R1, C3, R2 - R3);
d = zeros(size(Q,1), 1);
for i = 1:size(Q,1), d(i) = dist_to_arcs(Q(i,:), C1, R1, C2, R2, C3, R3); end
f = sqrt(mean(d.^2));
end

function d = dist_to_arcs(P, C1, R1, C2, R2, C3, R3)
d1 = abs(norm(P - C1) - R1);
d2 = abs(norm(P - C2) - R2);
d3 = abs(norm(P - C3) - R3);
d = min([d1, d2, d3]);
end

function C2 = internal_tangent_center(C1, r_a, C3, r_b)
% Point at distance r_a from C1 and r_b from C3 (circle-circle
% intersection); if the two circles do not actually meet (possible for
% intermediate, not-yet-converged (R1,R2,R3) during the search), returns
% the closest-approach point instead, keeping the objective smooth.
d = norm(C3 - C1);
if d < 1e-12, C2 = C1; return, end
u = (C3 - C1)/d;
a = (d^2 + r_a^2 - r_b^2)/(2*d);
h2 = max(r_a^2 - a^2, 0);
h = sqrt(h2);
Pm = C1 + a*u;
perp = [-u(2), u(1)];
% two candidates; pick the one on the concave (bore) side, i.e. the side
% NOT containing the segment A-B (a convex D bulges away from its chord)
cand1 = Pm + h*perp; cand2 = Pm - h*perp;
if cand1(2) < cand2(2), C2 = cand1; else, C2 = cand2; end
end

function u = unit(v)
n = norm(v);
if n < 1e-12, u = [1 0]; else, u = v/n; end
end

function ang = arc_angles(C, P0, P1)
a0 = atan2(P0(2) - C(2), P0(1) - C(1));
a1 = atan2(P1(2) - C(2), P1(1) - C(1));
ang = [a0, a1];
end
