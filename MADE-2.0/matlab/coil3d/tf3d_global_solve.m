function res = tf3d_global_solve(m)
%TF3D_GLOBAL_SOLVE Linear static solution of the 3D beam + shell global TF model.
%
%   res = TF3D_GLOBAL_SOLVE(m) solves the structural model built by
%   TF3D_GLOBAL_MODEL (the MATLAB counterpart of the ANSYS model STR_360):
%
%   - beams: 2-node 3D Timoshenko beams (exact cubic/shear-flexible
%     stiffness, 6 DOF per node), local x from node I to J, local z in the
%     plane of x and the orientation vector (as BEAM188 with node K), local
%     y = z x x. Section: A, Iy, Iz, J, Asy, Asz (m.sec, one per model).
%   - shells: 4-node flat shells (as SHELL181): membrane = bilinear quad
%     with Wilson incompatible modes (QM6, Taylor's correction, statically
%     condensed), plate = Mindlin with MITC4 transverse shear, drilling
%     rotation with a small penalty stiffness. Local x from node 1 to 2
%     (ANSYS element x axis I -> J), normal from the diagonals.
%   - supports: node, local frame (3 x 3, rows = axes in global), 1 x 6
%     mask of the fixed DOF in that frame (penalty-free: the constrained
%     node's DOF are rotated to the frame and removed).
%
%   m fields: X (N x 3), beam (E x 2), beam_k (E x 3 orientation vectors),
%   sec (struct E, nu, A, Iy, Iz, J, Asy, Asz), shell (S x 4), shell_t (S x
%   1), shell_E, shell_nu, sup (struct array: node, frame, mask), F (N x 6
%   nodal loads in global axes [Fx Fy Fz Mx My Mz]).
%
%   res: U (N x 6, global), beam_f (E x 12 end forces in local axes, forces
%   ON the element: [N_I Vy_I Vz_I T_I My_I Mz_I, N_J ... Mz_J] with
%   N_I > 0 = tension, as SMISC 1-6 / 14-19 of BEAM188 up to signs),
%   shell_r (S x 8 [N11 N22 N12 M11 M22 M12 Q13 Q23] at the centroid, local
%   axes, per unit length), reactions (struct per support: force/moment in
%   the support frame), info.

N = size(m.X, 1); nd = 6*N;
I = zeros(0,1); J = I; V = I;
nb = size(m.beam, 1);
Tb = cell(nb, 1); Kb = cell(nb, 1);
for e = 1:nb
    [kl, T] = beam_k(m.X(m.beam(e,1),:), m.X(m.beam(e,2),:), m.beam_k(e,:), m.sec);
    ke = T'*kl*T;
    d = [6*m.beam(e,1)-5:6*m.beam(e,1), 6*m.beam(e,2)-5:6*m.beam(e,2)];
    [ii, jj] = ndgrid(d, d);
    I = [I; ii(:)]; J = [J; jj(:)]; V = [V; ke(:)]; %#ok<AGROW>
    Tb{e} = T; Kb{e} = kl;
end
ns = size(m.shell, 1);
SH = cell(ns, 1);
for e = 1:ns
    s = shell_k(m.X(m.shell(e,:),:), m.shell_t(e), m.shell_E, m.shell_nu);
    d = reshape([6*m.shell(e,:)-5; 6*m.shell(e,:)-4; 6*m.shell(e,:)-3; ...
        6*m.shell(e,:)-2; 6*m.shell(e,:)-1; 6*m.shell(e,:)], 1, []);
    [ii, jj] = ndgrid(d, d);
    I = [I; ii(:)]; J = [J; jj(:)]; V = [V; s.K(:)]; %#ok<AGROW>
    s.dofs = d; SH{e} = s;
end
K = sparse(I, J, V, nd, nd);
K = (K + K')/2;
f = reshape(m.F', [], 1);
% supports: rotate the DOF of the supported nodes to their frame
R = speye(nd);
fixed = false(nd, 1);
for q = 1:numel(m.sup)
    n = m.sup(q).node; Fr = m.sup(q).frame;
    d = 6*n-5:6*n;
    R(d, d) = blkdiag(Fr', Fr');            % u_global = R * u_frame
    fixed(d(logical(m.sup(q).mask))) = true;
end
Kr = R'*K*R; fr = R'*f;
free = ~fixed;
ur = zeros(nd, 1);
Kff = Kr(free, free);
[L, flag, P] = chol(Kff, 'vector');
if flag == 0
    x = zeros(nnz(free), 1);
    b = fr(free);
    x(P) = L\(L'\b(P));
    ur(free) = x;
else
    ur(free) = Kff\fr(free);
end
u = R*ur;
res.U = reshape(u, 6, N)';
rr = Kr*ur - fr;
res.reactions = struct('node', {}, 'frame_force', {});
for q = 1:numel(m.sup)
    d = 6*m.sup(q).node-5:6*m.sup(q).node;
    res.reactions(q).node = m.sup(q).node;
    res.reactions(q).frame_force = rr(d)';
end
res.beam_f = zeros(nb, 12);
for e = 1:nb
    d = [6*m.beam(e,1)-5:6*m.beam(e,1), 6*m.beam(e,2)-5:6*m.beam(e,2)];
    fl = Kb{e}*(Tb{e}*u(d));
    res.beam_f(e,:) = [-fl(1:6)' fl(7:12)'];
    res.beam_f(e,1) = -fl(1);                % N_I > 0 tension
end
res.shell_r = zeros(ns, 8);
for e = 1:ns
    res.shell_r(e,:) = shell_resultants(SH{e}, u(SH{e}.dofs));
end
res.info = struct('ndof', nd, 'n_free', nnz(free), 'residual', norm(Kff*ur(free) - fr(free))/max(norm(fr(free)), eps));
end

% ------------------------------------------------------------------------
function [kl, T] = beam_k(xi, xj, kvec, s)
d = xj - xi; L = norm(d); ex = d/L;
ez = kvec - dot(kvec, ex)*ex; ez = ez/norm(ez);
ey = cross(ez, ex);
lam = [ex; ey; ez];
T = blkdiag(lam, lam, lam, lam);
E = s.E; G = E/(2*(1 + s.nu));
py = 12*E*s.Iz/(G*s.Asy*L^2);           % bending in x-y (about z), shear in y
pz = 12*E*s.Iy/(G*s.Asz*L^2);           % bending in x-z (about y), shear in z
kl = zeros(12);
a = E*s.A/L; t = G*s.J/L;
kl([1 7],[1 7]) = a*[1 -1; -1 1];
kl([4 10],[4 10]) = t*[1 -1; -1 1];
% x-y plane: v (2,8), theta_z (6,12)
c = E*s.Iz/(L^3*(1+py));
B = c*[12 6*L -12 6*L; 6*L (4+py)*L^2 -6*L (2-py)*L^2; -12 -6*L 12 -6*L; 6*L (2-py)*L^2 -6*L (4+py)*L^2];
kl([2 6 8 12],[2 6 8 12]) = B;
% x-z plane: w (3,9), theta_y (5,11) (theta_y = -dw/dx)
c = E*s.Iy/(L^3*(1+pz));
B = c*[12 -6*L -12 -6*L; -6*L (4+pz)*L^2 6*L (2-pz)*L^2; -12 6*L 12 6*L; -6*L (2-pz)*L^2 6*L (4+pz)*L^2];
kl([3 5 9 11],[3 5 9 11]) = B;
end

% ------------------------------------------------------------------------
function s = shell_k(Xg, t, E, nu)
% local frame
e1 = Xg(2,:) - Xg(1,:); e1 = e1/norm(e1);
n = cross(Xg(3,:) - Xg(1,:), Xg(4,:) - Xg(2,:)); n = n/norm(n);
e1 = e1 - dot(e1, n)*n; e1 = e1/norm(e1);
e2 = cross(n, e1);
lam = [e1; e2; n];
c0 = mean(Xg, 1);
xy = (Xg - c0)*[e1' e2'];
Dm = E*t/(1 - nu^2)*[1 nu 0; nu 1 0; 0 0 (1-nu)/2];
Db = E*t^3/(12*(1 - nu^2))*[1 nu 0; nu 1 0; 0 0 (1-nu)/2];
Ds = 5/6*E*t/(2*(1 + nu))*eye(2);
gp = [-1 1]/sqrt(3);
% membrane QM6: u = [u1 v1 ... u4 v4], incompatible a (4)
Kuu = zeros(8); Kua = zeros(8,4); Kaa = zeros(4);
[J0, ~, dJ0] = jac(xy, 0, 0);
Kb = zeros(12);                     % plate: [w th_x th_y] x 4
for a = gp
    for b = gp
        [Jm, dN, detJ] = jac(xy, a, b);
        Bm = zeros(3, 8);
        Bm(1, 1:2:8) = dN(1,:); Bm(2, 2:2:8) = dN(2,:);
        Bm(3, 1:2:8) = dN(2,:); Bm(3, 2:2:8) = dN(1,:);
        % incompatible modes (1 - xi^2), (1 - eta^2), Jacobian at centre
        dP = J0\[-2*a 0; 0 -2*b];       % d/dx, d/dy of the two modes
        dP = dP*dJ0/detJ;                % Taylor correction: integral of B_a vanishes (patch test)
        Ba = zeros(3, 4);
        Ba(1, [1 2]) = dP(1,:); Ba(2, [3 4]) = dP(2,:);
        Ba(3, [1 2]) = dP(2,:); Ba(3, [3 4]) = dP(1,:);
        Kuu = Kuu + Bm'*Dm*Bm*detJ; Kua = Kua + Bm'*Dm*Ba*detJ; Kaa = Kaa + Ba'*Dm*Ba*detJ;
        % plate bending: beta_x = th_y, beta_y = -th_x
        Bb = zeros(3, 12);
        Bb(1, 3:3:12) = dN(1,:);
        Bb(2, 2:3:12) = -dN(2,:);
        Bb(3, 3:3:12) = dN(2,:); Bb(3, 2:3:12) = -dN(1,:);
        Kb = Kb + Bb'*Db*Bb*detJ;
        % MITC4 transverse shear
        Bs = mitc4_Bs(xy, a, b, Jm);
        Kb = Kb + Bs'*Ds*Bs*detJ;
    end
end
Km = Kuu - Kua*(Kaa\Kua');
A = polyarea(xy(:,1), xy(:,2));
kdr = 1e-3*E*t*A/4;                     % drilling penalty
Kl = zeros(24);
for i = 1:4
    for j = 1:4
        Kl(6*i-5:6*i-4, 6*j-5:6*j-4) = Km(2*i-1:2*i, 2*j-1:2*j);
        Kl([6*i-3 6*i-2 6*i-1], [6*j-3 6*j-2 6*j-1]) = Kb(3*i-2:3*i, 3*j-2:3*j);
    end
    Kl(6*i, 6*i) = kdr;
end
T = kron(eye(8), lam);
s.K = T'*Kl*T;
s.T = T; s.xy = xy; s.Dm = Dm; s.Db = Db; s.Ds = Ds; s.Kua = Kua; s.Kaa = Kaa;
s.J0 = J0;
end

function [Jm, dN, detJ] = jac(xy, a, b)
dNn = 0.25*[-(1-b) (1-b) (1+b) -(1+b); -(1-a) -(1+a) (1+a) (1-a)];
Jm = dNn*xy;
detJ = det(Jm);
dN = Jm\dNn;
end

function Bs = mitc4_Bs(xy, a, b, Jm)
% covariant shear strains tied at the edge mid-points (Bathe-Dvorkin)
pts = [0 -1; 0 1; -1 0; 1 0];          % A(eta=-1), C(eta=+1) for e_xi; B(xi=-1), D(xi=+1) for e_eta
G = zeros(4, 12);
for q = 1:4
    xq = pts(q,1); yq = pts(q,2);
    Nq = 0.25*[(1-xq)*(1-yq) (1+xq)*(1-yq) (1+xq)*(1+yq) (1-xq)*(1+yq)];
    dNn = 0.25*[-(1-yq) (1-yq) (1+yq) -(1+yq); -(1-xq) -(1+xq) (1+xq) (1-xq)];
    Jq = dNn*xy;
    % gamma = dw/dx + beta, beta_x = th_y, beta_y = -th_x
    % covariant: gamma_xi = w,xi + beta . x,xi
    if q <= 2, r = 1; else, r = 2; end
    row = zeros(1, 12);
    row(1:3:12) = dNn(r,:);
    row(3:3:12) = Nq*Jq(r,1);
    row(2:3:12) = -Nq*Jq(r,2);
    G(q,:) = row;
end
gxi  = 0.5*(1-b)*G(1,:) + 0.5*(1+b)*G(2,:);
geta = 0.5*(1-a)*G(3,:) + 0.5*(1+a)*G(4,:);
Bs = Jm\[gxi; geta];
end

function r = shell_resultants(s, ug)
ul = s.T*ug;
um = zeros(8,1); up = zeros(12,1);
for i = 1:4
    um(2*i-1:2*i) = ul(6*i-5:6*i-4);
    up(3*i-2:3*i) = ul([6*i-3 6*i-2 6*i-1]);
end
[Jm, dN] = jac(s.xy, 0, 0);
Bm = zeros(3, 8);
Bm(1, 1:2:8) = dN(1,:); Bm(2, 2:2:8) = dN(2,:);
Bm(3, 1:2:8) = dN(2,:); Bm(3, 2:2:8) = dN(1,:);
Nr = s.Dm*Bm*um;                       % incompatible modes vanish at the centre
Bb = zeros(3, 12);
Bb(1, 3:3:12) = dN(1,:); Bb(2, 2:3:12) = -dN(2,:);
Bb(3, 3:3:12) = dN(2,:); Bb(3, 2:3:12) = -dN(1,:);
Mr = s.Db*Bb*up;
Qr = s.Ds*mitc4_Bs(s.xy, 0, 0, Jm)*up;
r = [Nr' Mr' Qr'];
end
