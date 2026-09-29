function out = tf3d_centreline_loads(P, WPH, WPW, NI, n_TF, pf, opts)
%TF3D_CENTRELINE_LOADS Field and nodal EM forces along a TF coil centreline.
%
%   out = TF3D_CENTRELINE_LOADS(P, WPH, WPW, NI, n_TF, pf, opts) computes
%   the magnetic field at the points P of the centreline of TF coil 1 and
%   the nodal Lorentz forces on them, with the same model as the ANSYS 3D
%   magnetic model (MAG_360: SOURC36 bars/arcs of the WP section along the
%   centreline, all n_TF coils, PF/CS/plasma rings, BIOT):
%
%   - P: M x 3 points of coil 1's centreline in order along the loop
%     (closed: the last point connects back to the first), coil 1 in the
%     plane x = 0 with y the major radius and z vertical (ANSYS CSYS 0);
%     the coils are at phi = 90 deg + k*360/n_TF.
%   - WPH, WPW: in-plane (radial) height and toroidal width of the
%     rectangular winding-pack section [m], scalars or one value per point
%     (the section of the source segment starting at that point); NI:
%     total current [A] (signed: positive flows from point k to point k+1).
%   - pf (optional): K x 5 [rc zc dr dz I] PF/CS/plasma rings coaxial with
%     the machine axis, rectangular section dr x dz, current I [A].
%   - opts: n_gauss (section Gauss points per direction, default 4),
%     n_ring (segments per PF ring, default 360), point_block/source_block,
%     eval_points (K x 3, optional: field also at these points, out.B_eval).
%
%   The rectangular section is represented by n_gauss x n_gauss filaments
%   that follow the centreline at the Gauss offsets (in-plane normal x
%   toroidal), each carrying NI/n_gauss^2 times the product of the Gauss
%   weights: the field of a uniformly distributed current, finite at the
%   centreline (no filament passes through it).
%
%   out: B (M x 3) [T], Bmag, F (M x 3) nodal forces [N] computed as the
%   APDL post-processing of MAG_360 (F = NI * dl x B with dl from central
%   differences, one-sided at the two ends), s (curvilinear abscissa).

if nargin < 6, pf = zeros(0, 5); end
if nargin < 7, opts = struct(); end
ng = getf(opts, 'n_gauss', 4); nr = getf(opts, 'n_ring', 360);
bs = struct('point_block', getf(opts, 'point_block', 128), 'source_block', getf(opts, 'source_block', 512));
M = size(P, 1);
% in-plane unit normal of the centreline (coil 1 plane x = 0) and the
% toroidal direction (x)
Pn = [P(end,:); P; P(1,:)];
t = Pn(3:end,:) - Pn(1:end-2,:); t = t./sqrt(sum(t.^2, 2));
n = [zeros(M,1), -t(:,3), t(:,2)];             % in-plane normal (rotated tangent)
e_tor = repmat([1 0 0], M, 1);
[xg, wg] = gauss_legendre(ng);
A = zeros(0,3); C = zeros(0,3); I = zeros(0,1);
for a = 1:ng
    for b = 1:ng
        Q = P + (xg(a)*WPH(:)/2).*n + (xg(b)*WPW(:)/2).*e_tor;
        Q2 = [Q(2:end,:); Q(1,:)];
        A = [A; Q]; C = [C; Q2]; %#ok<AGROW>
        I = [I; repmat(NI*wg(a)*wg(b)/4, M, 1)]; %#ok<AGROW>
    end
end
% other coils: rotation about z by k*2*pi/n_TF
As = A; Cs = C; Is = I;
for k = 1:n_TF-1
    th = k*2*pi/n_TF; R = [cos(th) -sin(th) 0; sin(th) cos(th) 0; 0 0 1];
    As = [As; A*R']; Cs = [Cs; C*R']; Is = [Is; I]; %#ok<AGROW>
end
% PF / CS / plasma rings
for k = 1:size(pf, 1)
    [xr, wr] = gauss_legendre(2);
    for a = 1:2
        for b = 1:2
            r = pf(k,1) + xr(a)*pf(k,3)/2; z = pf(k,2) + xr(b)*pf(k,4)/2;
            ph = (0:nr)'*2*pi/nr;
            Q = [r*cos(ph), r*sin(ph), z*ones(nr+1,1)];
            As = [As; Q(1:end-1,:)]; Cs = [Cs; Q(2:end,:)]; %#ok<AGROW>
            Is = [Is; repmat(pf(k,5)*wr(a)*wr(b)/4, nr, 1)]; %#ok<AGROW>
        end
    end
end
B = biot_savart_segments(P, As, Cs, Is, bs);
% nodal forces as the MAG_360 post-processing
dl = zeros(M, 3);
dl(1,:) = P(2,:) - P(1,:);
dl(2:M-1,:) = (P(3:M,:) - P(1:M-2,:))/2;
dl(M,:) = P(M,:) - P(M-1,:);
F = NI*cross(dl, B, 2);
out = struct('B', B, 'Bmag', sqrt(sum(B.^2, 2)), 'F', F, ...
    's', [0; cumsum(sqrt(sum(diff(P).^2, 2)))], 'n_sources', size(As, 1));
if isfield(opts, 'eval_points') && ~isempty(opts.eval_points)
    out.B_eval = biot_savart_segments(opts.eval_points, As, Cs, Is, bs);
end
end

function v = getf(s, f, d)
if isfield(s, f), v = s.(f); else, v = d; end
end

function [x, w] = gauss_legendre(n)
% Gauss-Legendre points and weights on [-1, 1]
b = (1:n-1)./sqrt(4*(1:n-1).^2 - 1);
[V, D] = eig(diag(b, 1) + diag(b, -1));
[x, i] = sort(diag(D)); w = 2*V(1,i)'.^2;
end
